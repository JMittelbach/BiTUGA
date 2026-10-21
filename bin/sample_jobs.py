"""Bounded sample scheduling and cancellable matcher/postprocessor pipelines."""
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
import subprocess
import tempfile
import threading


def run_sample_jobs(sample_ids, worker, jobs):
    """Call worker(sample_id, stop_event); stop scheduling on any failure."""
    stop = threading.Event()
    remaining = iter(sample_ids)
    executor = ThreadPoolExecutor(max_workers=jobs)
    pending = set()
    try:
        for _ in range(jobs):
            sid = next(remaining, None)
            if sid is not None:
                pending.add(executor.submit(worker, sid, stop))
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            # Check every completed job before admitting replacements.
            for future in done:
                future.result()
            for _ in done:
                sid = next(remaining, None)
                if sid is not None:
                    pending.add(executor.submit(worker, sid, stop))
    finally:
        stop.set()
        for future in pending:
            future.cancel()
        executor.shutdown(wait=True, cancel_futures=True)


def run_pipeline(matcher_cmd, postprocess_cmd, cwd, stop):
    """Drain diagnostics to disk, cancel peers on failure, and reap children."""
    processes = []
    with tempfile.TemporaryFile() as matcher_err, tempfile.TemporaryFile() as post_err:
        try:
            if stop.is_set():
                raise RuntimeError('Sample processing cancelled')
            matcher = subprocess.Popen(matcher_cmd, stdout=subprocess.PIPE,
                                       stderr=matcher_err, cwd=cwd)
            processes.append(matcher)
            try:
                post = subprocess.Popen(postprocess_cmd, stdin=matcher.stdout,
                                        stderr=post_err, cwd=cwd)
                processes.append(post)
            finally:
                matcher.stdout.close()
            while True:
                for process, label, errors in (
                    (matcher, 'nt_mini_matcher', matcher_err),
                    (post, 'postprocess_minimatcher', post_err),
                ):
                    rc = process.poll()
                    if rc is not None and rc != 0:
                        errors.seek(0, 2)
                        errors.seek(max(0, errors.tell() - 16384))
                        detail = errors.read().decode(errors='replace').strip()
                        raise RuntimeError(f'{label} failed with exit code {rc}: {detail}')
                if matcher.returncode == 0 and post.returncode == 0:
                    return
                if stop.wait(0.1):
                    raise RuntimeError('Sample processing cancelled')
        finally:
            for process in processes:
                if process.poll() is None:
                    process.terminate()
            for process in processes:
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
