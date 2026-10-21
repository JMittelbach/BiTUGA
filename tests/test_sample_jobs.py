from pathlib import Path
import random
import subprocess
import sys
import tempfile
import threading
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
from sample_jobs import run_pipeline, run_sample_jobs


class SchedulingTests(unittest.TestCase):
    def test_jobs_overlap_without_exceeding_limit(self):
        barrier = threading.Barrier(2, timeout=3)
        lock = threading.Lock()
        active, peak = 0, 0
        completed = []

        def worker(sid, stop):
            nonlocal active, peak
            with lock:
                active += 1
                peak = max(peak, active)
            barrier.wait()
            with lock:
                active -= 1
                completed.append(sid)

        run_sample_jobs(range(6), worker, 2)
        self.assertEqual(peak, 2)
        self.assertEqual(sorted(completed), list(range(6)))

    def test_failure_cancels_active_peer_and_admits_no_more_samples(self):
        barrier = threading.Barrier(2, timeout=3)
        seen = []

        def worker(sid, stop):
            seen.append(sid)
            barrier.wait()
            if sid == 'bad':
                raise RuntimeError('failed sample')
            self.assertTrue(stop.wait(3))

        with self.assertRaisesRegex(RuntimeError, 'failed sample'):
            run_sample_jobs(['bad', 'peer', 'never'], worker, 2)
        self.assertEqual(set(seen), {'bad', 'peer'})

    def test_large_stderr_and_failure_do_not_deadlock(self):
        with tempfile.TemporaryDirectory() as work:
            matcher = [sys.executable, '-c',
                       "import sys; sys.stderr.write('x'*200000); sys.exit(7)"]
            post = [sys.executable, '-c',
                    'import sys,time; sys.stdin.read(); time.sleep(30)']
            started = time.monotonic()
            with self.assertRaisesRegex(RuntimeError, 'exit code 7'):
                run_pipeline(matcher, post, work, threading.Event())
            self.assertLess(time.monotonic() - started, 10)

    def test_postprocess_failure_terminates_producer(self):
        with tempfile.TemporaryDirectory() as work:
            matcher = [sys.executable, '-c', 'import time; time.sleep(30)']
            post = [sys.executable, '-c', 'import sys; sys.exit(9)']
            started = time.monotonic()
            with self.assertRaisesRegex(RuntimeError, 'exit code 9'):
                run_pipeline(matcher, post, work, threading.Event())
            self.assertLess(time.monotonic() - started, 10)


@unittest.skipUnless(all((ROOT / 'bin' / name).is_file() for name in
                        ('nt_mini_matcher.x', 'postprocess_minimatcher', 'merge_matches')),
                     'Build native tools first')
class NativeSampleJobsTests(unittest.TestCase):
    def test_serial_and_parallel_presence_outputs_are_identical(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            rng = random.Random(342)
            sequences = [''.join(rng.choices('ACGT', k=80)) for _ in range(2)]
            metadata = work / 'metadata.tsv'
            rows = ['path\tsample_id\ttrait\n']
            for sid, which in [('a', [0]), ('b', [1]), ('c', [0, 1])]:
                lanes = 2 if sid == 'c' else 1
                for lane in range(lanes):
                    for mate in (1, 2):
                        fastq = work / f'{sid}_L{lane}_R{mate}.fq'
                        fastq.write_text(''.join(
                            f'@{sid}_{lane}_{idx}_{rep}/{mate}\n{sequences[idx]}\n+\n' + 'I'*80 + '\n'
                            for idx in which for rep in range(4 // lanes)))
                        rows.append(f'{fastq}\t{sid}\t{int(sid == "b")}\n')
            metadata.write_text(''.join(rows))
            outputs = []
            counts = []
            for jobs in (1, 2):
                out = work / f'jobs{jobs}'
                stats = out / 'results/statistics'
                stats.mkdir(parents=True)
                (out / 'results/unitigs.fa').write_text(''.join(
                    f'>{idx}\n{seq}\n' for idx, seq in enumerate(sequences)))
                (stats / 'unitig_statistics.txt').write_text(
                    'k\t31\nn_unitigs\t2\nmedian_len\t80\nmax_len\t80\n'
                    '# length_histogram\nlength\tcount\n80\t2\n')
                result = subprocess.run(
                    [sys.executable, '-B', str(ROOT / 'bin/orchestrator_4.py'),
                     '--outdir', str(out), '--metadata', str(metadata),
                     '--threads', '4', '--sample-jobs', str(jobs),
                     '--query-split-size', '4', '--min-mem-length', '32'],
                    capture_output=True, text=True, timeout=60)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                for sid in ('a', 'b', 'c'):
                    self.assertIn(f'[{sid}] DONE', result.stdout)
                counts.append([line.split('] [INFO] ')[1] for line in result.stdout.splitlines()
                               if 'PRESENT_UNITIGS\t' in line])
                outputs.append(out)
            self.assertEqual(sorted(counts[0]), sorted(counts[1]))
            self.assertEqual(sorted(counts[0]), ['[a] PRESENT_UNITIGS\t1',
                                               '[b] PRESENT_UNITIGS\t1',
                                               '[c] PRESENT_UNITIGS\t2'])
            self.assertEqual((outputs[0] / 'results/unitig_presence.tsv').read_bytes(),
                             (outputs[1] / 'results/unitig_presence.tsv').read_bytes())
