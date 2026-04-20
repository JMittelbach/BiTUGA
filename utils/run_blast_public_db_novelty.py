#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen qstart qend "
    "sstart send evalue bitscore qcovs staxids sscinames stitle"
)


@dataclass
class DatasetConfig:
    key: str
    label: str
    fasta: Path
    unmapped_ids: Path


def now_ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str, level: str = "INFO") -> None:
    print(f"[{now_ts()}] [{level}] {msg}", file=sys.stderr, flush=True)


def fmt_duration(seconds: float) -> str:
    sec = int(round(max(0.0, seconds)))
    h, rem = divmod(sec, 3600)
    m, s = divmod(rem, 60)
    if h > 0:
        return f"{h:d}h {m:02d}m {s:02d}s"
    return f"{m:d}m {s:02d}s"


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    repo_root = here.parent
    revision_dir = repo_root / "doc/paper/G3/revision"
    paper_figshare_dir = repo_root / "paper_figshare"

    ap = argparse.ArgumentParser(
        description=(
            "Run batch BLASTn against public DBs for Populus/Ginkgo trait-associated unitigs "
            "and summarize hit-vs-novel statistics, including unmapped-reference subsets."
        )
    )
    ap.add_argument(
        "--populus-fasta",
        type=Path,
        default=repo_root / "paper_figshare/populus/trait_associated_unitigs_populus.fa",
    )
    ap.add_argument(
        "--ginkgo-fasta",
        type=Path,
        default=repo_root / "paper_figshare/ginkgo/trait_associated_unitigs_ginkgo.fa",
    )
    ap.add_argument(
        "--populus-unmapped-ids",
        type=Path,
        default=repo_root / "paper_figshare/populus/statistics/populus_unmapped_unitig_ids_minimap.txt",
    )
    ap.add_argument(
        "--ginkgo-unmapped-ids",
        type=Path,
        default=repo_root / "paper_figshare/ginkgo/statistics/ginkgo_unmapped_unitig_ids_minimap.txt",
    )
    ap.add_argument("--outdir", type=Path, default=revision_dir / "blast_public_db_results")
    ap.add_argument(
        "--summary-txt",
        type=Path,
        default=revision_dir / "blast_public_db_novelty_summary.txt",
    )
    ap.add_argument(
        "--supplementary5-tsv",
        type=Path,
        default=paper_figshare_dir / "Supplementary_File_5.tsv",
        help="Output path for the human-readable Supplementary File 5 report.",
    )

    ap.add_argument("--blastn-bin", default="blastn")
    ap.add_argument("--db", default="nt")
    ap.add_argument("--task", default="megablast")
    ap.add_argument("--max-target-seqs", type=int, default=20)
    ap.add_argument("--max-hsps", type=int, default=1)
    ap.add_argument("--batch-size", type=int, default=50)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--no-remote", action="store_true", help="Disable BLAST -remote mode.")
    ap.add_argument(
        "--parse-only",
        action="store_true",
        help="Skip BLAST execution and only parse existing combined/batch outputs.",
    )
    ap.add_argument(
        "--resume",
        action="store_true",
        help="Legacy alias; resume is now enabled by default unless --no-resume is set.",
    )
    ap.add_argument(
        "--no-resume",
        action="store_true",
        help="Disable internal auto-resume and rerun all batches from scratch.",
    )
    ap.add_argument(
        "--combined-hits-name",
        default="{dataset}_blast_hits_all_batches.tsv",
        help="Output filename pattern for cumulative BLAST rows per dataset.",
    )
    ap.add_argument(
        "--combined-significant-hits-name",
        default="{dataset}_blast_hits_significant_all_batches.tsv",
        help="Output filename pattern for cumulative significant BLAST rows per dataset.",
    )
    ap.add_argument(
        "--checkpoint-name",
        default="{dataset}_resume_checkpoint.json",
        help="Output filename pattern for per-dataset resume checkpoint JSON.",
    )
    ap.add_argument(
        "--cleanup-batch-tsv-after-combine",
        action="store_true",
        help=(
            "Deprecated no-op: per-batch BLAST TSV files are now deleted immediately "
            "after successful rolling accumulation."
        ),
    )
    ap.add_argument(
        "--cleanup-batch-fasta-after-combine",
        action="store_true",
        help=(
            "Deprecated no-op: per-batch FASTA files are now deleted immediately "
            "after successful rolling accumulation."
        ),
    )
    ap.add_argument("--max-retries", type=int, default=3)
    ap.add_argument("--retry-wait-sec", type=float, default=25.0)
    ap.add_argument("--delay-sec", type=float, default=3.0, help="Sleep between successful batch submissions.")
    ap.add_argument(
        "--limit-queries-per-dataset",
        type=int,
        default=0,
        help="For testing only: limit the number of query sequences processed per dataset.",
    )
    ap.add_argument(
        "--datasets",
        nargs="+",
        choices=("populus", "ginkgo"),
        default=("populus", "ginkgo"),
        help="Datasets to process (default: both).",
    )

    ap.add_argument("--hit-evalue", type=float, default=1e-10)
    ap.add_argument("--hit-qcovs", type=float, default=70.0)
    ap.add_argument("--hit-pident", type=float, default=85.0)
    ap.add_argument(
        "--hit-length",
        type=int,
        default=0,
        help="Absolute minimum alignment length (bp). Use 0 to disable absolute cutoff.",
    )
    ap.add_argument(
        "--hit-min-aln-frac",
        type=float,
        default=0.90,
        help="Minimum aligned fraction of query length (0..1).",
    )

    args = ap.parse_args()
    if args.batch_size <= 0:
        ap.error("--batch-size must be > 0")
    if args.max_target_seqs <= 0:
        ap.error("--max-target-seqs must be > 0")
    if args.max_hsps <= 0:
        ap.error("--max-hsps must be > 0")
    if args.threads <= 0:
        ap.error("--threads must be > 0")
    if args.max_retries <= 0:
        ap.error("--max-retries must be > 0")
    if args.hit_length < 0:
        ap.error("--hit-length must be >= 0")
    if not (0.0 <= args.hit_min_aln_frac <= 1.0):
        ap.error("--hit-min-aln-frac must be between 0 and 1")
    if args.resume and args.no_resume:
        ap.error("--resume and --no-resume are mutually exclusive")
    args.auto_resume = bool(args.resume or not args.no_resume)
    args.repo_root = repo_root
    return args


def read_fasta_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    cur_header = None
    cur_seq: list[str] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_header is not None:
                    records.append((cur_header, "".join(cur_seq)))
                cur_header = line[1:].strip()
                cur_seq = []
            else:
                cur_seq.append(line)
    if cur_header is not None:
        records.append((cur_header, "".join(cur_seq)))
    return records


def fasta_id(header: str) -> str:
    return header.split()[0]


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def load_id_set(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open() as fh:
        for raw in fh:
            v = raw.strip()
            if v:
                ids.add(v)
    return ids


def run_blast_batch(
    args: argparse.Namespace,
    query_fasta: Path,
    out_tsv: Path,
) -> None:
    remote = not args.no_remote
    cmd = [
        args.blastn_bin,
        "-query",
        str(query_fasta),
        "-db",
        args.db,
        "-task",
        args.task,
        "-max_target_seqs",
        str(args.max_target_seqs),
        "-max_hsps",
        str(args.max_hsps),
        "-outfmt",
        OUTFMT,
        "-out",
        str(out_tsv),
    ]
    if remote:
        cmd.append("-remote")
    else:
        cmd.extend(["-num_threads", str(args.threads)])

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    tmp_out = out_tsv.with_suffix(out_tsv.suffix + ".tmp")

    last_err = None
    for attempt in range(1, args.max_retries + 1):
        if tmp_out.exists():
            tmp_out.unlink()
        run_cmd = cmd.copy()
        run_cmd[run_cmd.index(str(out_tsv))] = str(tmp_out)
        attempt_start = time.time()
        log(f"BLAST start batch={query_fasta.name} attempt={attempt}/{args.max_retries}")
        proc = subprocess.run(run_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if proc.returncode == 0 and tmp_out.exists():
            tmp_out.replace(out_tsv)
            log(
                f"BLAST done batch={query_fasta.name} "
                f"attempt={attempt}/{args.max_retries} elapsed={fmt_duration(time.time() - attempt_start)}"
            )
            return
        last_err = proc.stderr.strip() if proc.stderr else f"exit={proc.returncode}"
        log(f"BLAST failed batch={query_fasta.name}: {last_err}", level="WARN")
        if attempt < args.max_retries:
            time.sleep(args.retry_wait_sec)

    raise RuntimeError(f"BLAST failed for {query_fasta.name}: {last_err}")


def parse_blast_tsv(
    tsv_paths: list[Path],
    query_len_by_id: dict[str, int],
    hit_evalue: float,
    hit_qcovs: float,
    hit_pident: float,
    hit_length: int,
    hit_min_aln_frac: float,
) -> tuple[set[str], set[str], int]:
    any_hit_ids: set[str] = set()
    sig_hit_ids: set[str] = set()
    row_count = 0
    for path in sorted(tsv_paths):
        if not path.exists():
            continue
        with path.open() as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                row_count += 1
                cols = line.split("\t")
                if len(cols) < 13:
                    continue
                qseqid = cols[0]
                any_hit_ids.add(qseqid)
                if is_significant_hit_cols(
                    cols,
                    query_len_by_id=query_len_by_id,
                    hit_evalue=hit_evalue,
                    hit_qcovs=hit_qcovs,
                    hit_pident=hit_pident,
                    hit_length=hit_length,
                    hit_min_aln_frac=hit_min_aln_frac,
                ):
                    sig_hit_ids.add(qseqid)
    return any_hit_ids, sig_hit_ids, row_count


def is_significant_hit_cols(
    cols: list[str],
    query_len_by_id: dict[str, int],
    hit_evalue: float,
    hit_qcovs: float,
    hit_pident: float,
    hit_length: int,
    hit_min_aln_frac: float,
) -> bool:
    if len(cols) < 13:
        return False
    try:
        qseqid = cols[0]
        pident = float(cols[2])
        aln_len = int(cols[3])
        evalue = float(cols[10])
        qcovs = float(cols[12])
    except ValueError:
        return False
    min_len = max(0, hit_length)
    qlen = query_len_by_id.get(qseqid)
    if qlen is not None and hit_min_aln_frac > 0.0:
        min_len = max(min_len, int(math.ceil(hit_min_aln_frac * qlen)))
    return (
        evalue <= hit_evalue
        and qcovs >= hit_qcovs
        and pident >= hit_pident
        and aln_len >= min_len
    )


def write_combined_hits(
    tsv_paths: list[Path],
    out_all: Path,
    out_sig: Path,
    query_len_by_id: dict[str, int],
    hit_evalue: float,
    hit_qcovs: float,
    hit_pident: float,
    hit_length: int,
    hit_min_aln_frac: float,
) -> tuple[int, int]:
    out_all.parent.mkdir(parents=True, exist_ok=True)
    out_sig.parent.mkdir(parents=True, exist_ok=True)
    n_all = 0
    n_sig = 0
    with out_all.open("w") as fh_all, out_sig.open("w") as fh_sig:
        for path in sorted(tsv_paths):
            if not path.exists():
                continue
            with path.open() as fh:
                for raw in fh:
                    line = raw.strip()
                    if not line:
                        continue
                    n_all += 1
                    fh_all.write(line + "\n")
                    cols = line.split("\t")
                    if is_significant_hit_cols(
                        cols,
                        query_len_by_id=query_len_by_id,
                        hit_evalue=hit_evalue,
                        hit_qcovs=hit_qcovs,
                        hit_pident=hit_pident,
                        hit_length=hit_length,
                        hit_min_aln_frac=hit_min_aln_frac,
                    ):
                        n_sig += 1
                        fh_sig.write(line + "\n")
    return n_all, n_sig


def parse_combined_hits(
    out_all: Path,
    out_sig: Path,
) -> tuple[set[str], set[str], int, int]:
    any_hit_ids: set[str] = set()
    sig_hit_ids: set[str] = set()
    n_all = 0
    n_sig = 0

    if out_all.exists():
        with out_all.open() as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                n_all += 1
                cols = line.split("\t")
                if cols:
                    any_hit_ids.add(cols[0])

    if out_sig.exists():
        with out_sig.open() as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                n_sig += 1
                cols = line.split("\t")
                if cols:
                    sig_hit_ids.add(cols[0])

    return any_hit_ids, sig_hit_ids, n_all, n_sig


def append_batch_to_combined(
    batch_tsv: Path,
    out_all: Path,
    out_sig: Path,
    query_len_by_id: dict[str, int],
    hit_evalue: float,
    hit_qcovs: float,
    hit_pident: float,
    hit_length: int,
    hit_min_aln_frac: float,
) -> tuple[set[str], set[str], int, int]:
    out_all.parent.mkdir(parents=True, exist_ok=True)
    out_sig.parent.mkdir(parents=True, exist_ok=True)

    batch_any_ids: set[str] = set()
    batch_sig_ids: set[str] = set()
    n_all = 0
    n_sig = 0

    with batch_tsv.open() as fh_in, out_all.open("a") as fh_all, out_sig.open("a") as fh_sig:
        for raw in fh_in:
            line = raw.strip()
            if not line:
                continue
            n_all += 1
            fh_all.write(line + "\n")
            cols = line.split("\t")
            if cols:
                batch_any_ids.add(cols[0])
            if is_significant_hit_cols(
                cols,
                query_len_by_id=query_len_by_id,
                hit_evalue=hit_evalue,
                hit_qcovs=hit_qcovs,
                hit_pident=hit_pident,
                hit_length=hit_length,
                hit_min_aln_frac=hit_min_aln_frac,
            ):
                n_sig += 1
                fh_sig.write(line + "\n")
                if cols:
                    batch_sig_ids.add(cols[0])

    return batch_any_ids, batch_sig_ids, n_all, n_sig


def write_id_list(path: Path, ids: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for uid in sorted(ids):
            fh.write(uid + "\n")


def cleanup_files(paths: Iterable[Path]) -> tuple[int, int]:
    removed = 0
    missing = 0
    for p in paths:
        if p.exists():
            p.unlink()
            removed += 1
        else:
            missing += 1
    return removed, missing


def pct(n: int, d: int) -> str:
    if d == 0:
        return "NA"
    return f"{(100.0 * n / d):.2f}%"


def pct_num(n: int, d: int) -> str:
    if d == 0:
        return "NA"
    return f"{(100.0 * n / d):.6f}"


def process_dataset(args: argparse.Namespace, ds: DatasetConfig) -> dict[str, object]:
    records = read_fasta_records(ds.fasta)
    if args.limit_queries_per_dataset > 0:
        records = records[: args.limit_queries_per_dataset]
    query_len_by_id = {fasta_id(h): len(seq) for h, seq in records}
    all_ids = [fasta_id(h) for h, _ in records]
    all_id_set = set(all_ids)
    total_n = len(all_ids)
    if len(all_id_set) != total_n:
        raise RuntimeError(f"Duplicate FASTA IDs detected for {ds.label}; cannot classify robustly.")

    unmapped_ids_full = load_id_set(ds.unmapped_ids)
    unmapped_ids = unmapped_ids_full & all_id_set
    mapped_ids = all_id_set - unmapped_ids

    ds_dir = args.outdir / ds.key
    batch_fasta_dir = ds_dir / "batch_fasta"
    batch_tsv_dir = ds_dir / "batch_blast_tsv"
    batch_log_tsv = ds_dir / f"{ds.key}_batch_progress.tsv"
    checkpoint_json = ds_dir / args.checkpoint_name.format(dataset=ds.key)
    combined_all = ds_dir / args.combined_hits_name.format(dataset=ds.key)
    combined_sig = ds_dir / args.combined_significant_hits_name.format(dataset=ds.key)
    ds_dir.mkdir(parents=True, exist_ok=True)
    batch_fasta_dir.mkdir(parents=True, exist_ok=True)
    batch_tsv_dir.mkdir(parents=True, exist_ok=True)

    if not batch_log_tsv.exists():
        batch_log_tsv.write_text(
            "timestamp_utc\tdataset\tbatch\tstatus\tbatch_fasta\tbatch_tsv\telapsed_sec\tnote\n"
        )

    num_batches = int(math.ceil(total_n / args.batch_size)) if total_n > 0 else 0

    def append_batch_log(
        batch_no: int,
        status: str,
        batch_fa: Path,
        batch_tsv: Path,
        elapsed_sec: float,
        note: str,
    ) -> None:
        safe_note = str(note).replace("\t", " ").replace("\n", " ").strip()
        with batch_log_tsv.open("a") as lfh:
            lfh.write(
                f"{now_ts()}\t{ds.key}\t{batch_no:04d}\t{status}\t{batch_fa}\t{batch_tsv}\t{elapsed_sec:.3f}\t{safe_note}\n"
            )

    def write_checkpoint(
        last_started_batch: int,
        last_completed_batch: int,
        status: str,
        note: str,
        combined_up_to_batch: int,
        combined_rows_all: int,
        combined_rows_sig: int,
    ) -> None:
        payload = {
            "dataset": ds.key,
            "fasta": str(ds.fasta),
            "unmapped_ids": str(ds.unmapped_ids),
            "hit_evalue": float(args.hit_evalue),
            "hit_qcovs": float(args.hit_qcovs),
            "hit_pident": float(args.hit_pident),
            "hit_length": int(args.hit_length),
            "hit_min_aln_frac": float(args.hit_min_aln_frac),
            "batch_size": args.batch_size,
            "num_batches": num_batches,
            "last_started_batch": int(last_started_batch),
            "last_completed_batch": int(last_completed_batch),
            "combined_up_to_batch": int(combined_up_to_batch),
            "combined_rows_all": int(combined_rows_all),
            "combined_rows_sig": int(combined_rows_sig),
            "status": status,
            "note": note,
            "updated_at": now_ts(),
        }
        checkpoint_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    resume_from_batch = 1
    last_started_batch = 0
    last_completed_batch = 0
    combined_up_to_batch = 0
    checkpoint_loaded = False
    any_hit_ids, sig_hit_ids, row_count, combined_sig_rows = parse_combined_hits(
        combined_all, combined_sig
    )

    if not args.auto_resume and not args.parse_only:
        for p in (combined_all, combined_sig, checkpoint_json):
            if p.exists():
                p.unlink()
        removed_tsv = 0
        removed_fa = 0
        for p in batch_tsv_dir.glob("batch_*.tsv*"):
            p.unlink()
            removed_tsv += 1
        for p in batch_fasta_dir.glob("batch_*.fa"):
            p.unlink()
            removed_fa += 1
        log(
            f"{ds.key}: --no-resume -> cleaned previous artifacts "
            f"(combined/checkpoint + batch_tsv={removed_tsv}, batch_fasta={removed_fa})"
        )
        any_hit_ids, sig_hit_ids, row_count, combined_sig_rows = set(), set(), 0, 0

    if args.auto_resume and checkpoint_json.exists():
        try:
            cp = json.loads(checkpoint_json.read_text())
        except Exception as exc:
            log(f"{ds.key}: checkpoint unreadable -> ignore ({exc})", level="WARN")
        else:
            compatible = (
                cp.get("dataset") == ds.key
                and int(cp.get("batch_size", -1)) == args.batch_size
                and int(cp.get("num_batches", -1)) == num_batches
                and str(cp.get("fasta", "")) == str(ds.fasta)
                and str(cp.get("unmapped_ids", "")) == str(ds.unmapped_ids)
                and float(cp.get("hit_evalue", -1.0)) == float(args.hit_evalue)
                and float(cp.get("hit_qcovs", -1.0)) == float(args.hit_qcovs)
                and float(cp.get("hit_pident", -1.0)) == float(args.hit_pident)
                and int(cp.get("hit_length", -1)) == int(args.hit_length)
                and float(cp.get("hit_min_aln_frac", -1.0)) == float(args.hit_min_aln_frac)
            )
            if compatible:
                checkpoint_loaded = True
                last_started_batch = int(cp.get("last_started_batch", 0))
                last_completed_batch = int(cp.get("last_completed_batch", 0))
                combined_up_to_batch = int(cp.get("combined_up_to_batch", last_completed_batch))
                resume_from_batch = max(1, last_completed_batch + 1)
                log(
                    f"{ds.key}: checkpoint resume target batch "
                    f"{resume_from_batch:04d}/{num_batches:04d} "
                    f"(last_completed={last_completed_batch:04d})"
                )
            else:
                log(
                    f"{ds.key}: checkpoint metadata mismatch -> ignore and infer from files",
                    level="WARN",
                )

    if args.auto_resume and not checkpoint_loaded and (combined_all.exists() or combined_sig.exists()):
        raise RuntimeError(
            f"{ds.key}: found existing combined TSVs but no valid checkpoint. "
            "Delete combined outputs for a fresh run or rerun with --no-resume."
        )

    if checkpoint_loaded:
        if combined_up_to_batch < last_completed_batch:
            log(
                f"{ds.key}: checkpoint inconsistent (combined_up_to < last_completed). "
                f"Using combined_up_to={combined_up_to_batch:04d}.",
                level="WARN",
            )
            last_completed_batch = combined_up_to_batch
            resume_from_batch = max(1, last_completed_batch + 1)

        if last_completed_batch > 0 and row_count == 0 and not args.parse_only:
            existing_batch_tsv = sorted(batch_tsv_dir.glob("batch_*.tsv"))
            if existing_batch_tsv:
                log(
                    f"{ds.key}: combined TSVs missing but found {len(existing_batch_tsv)} batch TSVs; "
                    "reconstructing combined files",
                    level="WARN",
                )
                if combined_all.exists():
                    combined_all.unlink()
                if combined_sig.exists():
                    combined_sig.unlink()
                any_hit_ids, sig_hit_ids, row_count, combined_sig_rows = set(), set(), 0, 0
                reconstructed_max_batch = 0
                for path in existing_batch_tsv:
                    name = path.stem
                    if name.startswith("batch_"):
                        try:
                            reconstructed_max_batch = max(
                                reconstructed_max_batch, int(name.split("_", 1)[1])
                            )
                        except ValueError:
                            pass
                    b_any, b_sig, b_rows, b_sig_rows = append_batch_to_combined(
                        path,
                        out_all=combined_all,
                        out_sig=combined_sig,
                        query_len_by_id=query_len_by_id,
                        hit_evalue=args.hit_evalue,
                        hit_qcovs=args.hit_qcovs,
                        hit_pident=args.hit_pident,
                        hit_length=args.hit_length,
                        hit_min_aln_frac=args.hit_min_aln_frac,
                    )
                    any_hit_ids.update(b_any)
                    sig_hit_ids.update(b_sig)
                    row_count += b_rows
                    combined_sig_rows += b_sig_rows
                if reconstructed_max_batch > 0 and reconstructed_max_batch < last_completed_batch:
                    log(
                        f"{ds.key}: reconstructed only up to batch {reconstructed_max_batch:04d}; "
                        "adjusting checkpoint resume target",
                        level="WARN",
                    )
                    last_completed_batch = reconstructed_max_batch
                    combined_up_to_batch = min(combined_up_to_batch, reconstructed_max_batch)
                    resume_from_batch = max(1, last_completed_batch + 1)
            else:
                raise RuntimeError(
                    f"{ds.key}: checkpoint says completed batches exist, but combined TSV is empty/missing "
                    "and no batch TSVs are available for reconstruction."
                )

        if (
            last_started_batch > last_completed_batch
            and 1 <= last_started_batch <= num_batches
        ):
            interrupted = last_started_batch
            interrupted_tsv = batch_tsv_dir / f"batch_{interrupted:04d}.tsv"
            interrupted_tmp = interrupted_tsv.with_suffix(interrupted_tsv.suffix + ".tmp")
            interrupted_fa = batch_fasta_dir / f"batch_{interrupted:04d}.fa"
            removed = []
            for p in (interrupted_tsv, interrupted_tmp, interrupted_fa):
                if p.exists():
                    p.unlink()
                    removed.append(str(p))
            note = (
                f"detected interrupted batch {interrupted:04d}; "
                f"deleted stale files and rerun from this batch"
            )
            if removed:
                note += f" ({', '.join(removed)})"
            log(f"{ds.key}: {note}", level="WARN")
            append_batch_log(
                interrupted,
                "resume_requeue",
                interrupted_fa,
                interrupted_tsv,
                0.0,
                note,
            )
            resume_from_batch = min(resume_from_batch, interrupted)

    if args.parse_only:
        if row_count == 0 and combined_sig_rows == 0:
            batch_tsv_paths = sorted(batch_tsv_dir.glob("batch_*.tsv"))
            if batch_tsv_paths:
                if combined_all.exists():
                    combined_all.unlink()
                if combined_sig.exists():
                    combined_sig.unlink()
                for path in batch_tsv_paths:
                    b_any, b_sig, b_rows, b_sig_rows = append_batch_to_combined(
                        path,
                        out_all=combined_all,
                        out_sig=combined_sig,
                        query_len_by_id=query_len_by_id,
                        hit_evalue=args.hit_evalue,
                        hit_qcovs=args.hit_qcovs,
                        hit_pident=args.hit_pident,
                        hit_length=args.hit_length,
                        hit_min_aln_frac=args.hit_min_aln_frac,
                    )
                    any_hit_ids.update(b_any)
                    sig_hit_ids.update(b_sig)
                    row_count += b_rows
                    combined_sig_rows += b_sig_rows
                log(
                    f"{ds.key}: parse-only reconstructed combined files from {len(batch_tsv_paths)} batch TSVs"
                )
        write_checkpoint(
            last_started_batch=last_started_batch,
            last_completed_batch=last_completed_batch,
            status="parse_only",
            note="no BLAST execution",
            combined_up_to_batch=combined_up_to_batch,
            combined_rows_all=row_count,
            combined_rows_sig=combined_sig_rows,
        )
    else:
        if resume_from_batch > num_batches:
            log(f"{ds.key}: all batches already completed; nothing to run")
        else:
            log(f"{ds.key}: rolling accumulate resume from batch {resume_from_batch:04d}/{num_batches:04d}")

        for batch_no in range(resume_from_batch, num_batches + 1):
            start = (batch_no - 1) * args.batch_size
            end = min(total_n, start + args.batch_size)
            batch_records = records[start:end]
            batch_fa = batch_fasta_dir / f"batch_{batch_no:04d}.fa"
            batch_tsv = batch_tsv_dir / f"batch_{batch_no:04d}.tsv"
            batch_tmp = batch_tsv.with_suffix(batch_tsv.suffix + ".tmp")

            stale_removed = []
            for p in (batch_tsv, batch_tmp, batch_fa):
                if p.exists():
                    p.unlink()
                    stale_removed.append(str(p))
            if stale_removed:
                log(
                    f"{ds.key}: removed stale artifacts before batch {batch_no:04d}: "
                    f"{', '.join(stale_removed)}",
                    level="WARN",
                )

            write_fasta(batch_fa, batch_records)
            last_started_batch = batch_no
            write_checkpoint(
                last_started_batch=batch_no,
                last_completed_batch=last_completed_batch,
                status="running",
                note=f"started batch {batch_no:04d}",
                combined_up_to_batch=combined_up_to_batch,
                combined_rows_all=row_count,
                combined_rows_sig=combined_sig_rows,
            )
            append_batch_log(
                batch_no, "started", batch_fa, batch_tsv, 0.0, f"n_queries={len(batch_records)}"
            )
            t0 = time.time()
            log(f"{ds.key}: batch {batch_no:04d}/{num_batches:04d} start n_queries={len(batch_records)}")

            try:
                run_blast_batch(args, batch_fa, batch_tsv)
            except Exception as exc:
                elapsed = time.time() - t0
                append_batch_log(batch_no, "failed", batch_fa, batch_tsv, elapsed, str(exc))
                write_checkpoint(
                    last_started_batch=batch_no,
                    last_completed_batch=last_completed_batch,
                    status="failed",
                    note=str(exc),
                    combined_up_to_batch=combined_up_to_batch,
                    combined_rows_all=row_count,
                    combined_rows_sig=combined_sig_rows,
                )
                raise

            append_t0 = time.time()
            b_any, b_sig, b_rows, b_sig_rows = append_batch_to_combined(
                batch_tsv,
                out_all=combined_all,
                out_sig=combined_sig,
                query_len_by_id=query_len_by_id,
                hit_evalue=args.hit_evalue,
                hit_qcovs=args.hit_qcovs,
                hit_pident=args.hit_pident,
                hit_length=args.hit_length,
                hit_min_aln_frac=args.hit_min_aln_frac,
            )
            append_elapsed = time.time() - append_t0
            any_hit_ids.update(b_any)
            sig_hit_ids.update(b_sig)
            row_count += b_rows
            combined_sig_rows += b_sig_rows

            last_completed_batch = batch_no
            combined_up_to_batch = batch_no
            elapsed = time.time() - t0
            append_batch_log(
                batch_no,
                "done",
                batch_fa,
                batch_tsv,
                elapsed,
                f"rows={b_rows}; sig_rows={b_sig_rows}; append_elapsed={append_elapsed:.3f}s",
            )
            write_checkpoint(
                last_started_batch=batch_no,
                last_completed_batch=last_completed_batch,
                status="done",
                note=f"rows={b_rows}; sig_rows={b_sig_rows}",
                combined_up_to_batch=combined_up_to_batch,
                combined_rows_all=row_count,
                combined_rows_sig=combined_sig_rows,
            )

            removed = []
            for p in (batch_tsv, batch_tmp, batch_fa):
                if p.exists():
                    p.unlink()
                    removed.append(str(p))
            append_batch_log(
                batch_no,
                "cleanup",
                batch_fa,
                batch_tsv,
                0.0,
                f"removed_files={len(removed)}",
            )
            log(
                f"{ds.key}: batch {batch_no:04d}/{num_batches:04d} "
                f"done elapsed={fmt_duration(elapsed)} rows={b_rows} sig_rows={b_sig_rows}"
            )
            time.sleep(args.delay_sec)

        write_checkpoint(
            last_started_batch=max(last_started_batch, num_batches),
            last_completed_batch=num_batches,
            status="complete",
            note="all batches processed",
            combined_up_to_batch=num_batches,
            combined_rows_all=row_count,
            combined_rows_sig=combined_sig_rows,
        )

    if args.parse_only:
        any_hit_ids, sig_hit_ids, row_count, combined_sig_rows = parse_combined_hits(
            combined_all, combined_sig
        )

    novel_ids = all_id_set - sig_hit_ids
    unmapped_hit_ids = unmapped_ids & sig_hit_ids
    unmapped_novel_ids = unmapped_ids - sig_hit_ids

    all_query_ids_file = ds_dir / f"{ds.key}_all_query_ids.txt"
    hit_ids_file = ds_dir / f"{ds.key}_hit_ids_significant.txt"
    novel_ids_file = ds_dir / f"{ds.key}_novel_ids_no_significant_hit.txt"
    unmapped_hit_ids_file = ds_dir / f"{ds.key}_unmapped_hit_ids_significant.txt"
    unmapped_novel_ids_file = ds_dir / f"{ds.key}_unmapped_novel_ids.txt"

    write_id_list(all_query_ids_file, all_id_set)
    write_id_list(hit_ids_file, sig_hit_ids)
    write_id_list(novel_ids_file, novel_ids)
    write_id_list(unmapped_hit_ids_file, unmapped_hit_ids)
    write_id_list(unmapped_novel_ids_file, unmapped_novel_ids)

    return {
        "dataset_key": ds.key,
        "dataset": ds.label,
        "input_fasta": ds.fasta,
        "input_unmapped_ids": ds.unmapped_ids,
        "total_unitigs": total_n,
        "mapped_unitigs_ref": len(mapped_ids),
        "unmapped_unitigs_ref": len(unmapped_ids),
        "queries_with_any_blast_row": len(any_hit_ids),
        "queries_with_significant_hit": len(sig_hit_ids),
        "queries_without_significant_hit_novel": len(novel_ids),
        "unmapped_with_significant_hit": len(unmapped_hit_ids),
        "unmapped_novel": len(unmapped_novel_ids),
        "blast_rows_total": row_count,
        "num_batches": num_batches,
        "batch_size": args.batch_size,
        "batch_progress_log": batch_log_tsv,
        "resume_checkpoint_json": checkpoint_json,
        "combined_hits_all_tsv": combined_all,
        "combined_hits_significant_tsv": combined_sig,
        "all_query_ids_file": all_query_ids_file,
        "hit_ids_file": hit_ids_file,
        "novel_ids_file": novel_ids_file,
        "unmapped_hit_ids_file": unmapped_hit_ids_file,
        "unmapped_novel_ids_file": unmapped_novel_ids_file,
        "result_dir": ds_dir,
    }


def write_summary_txt(path: Path, args: argparse.Namespace, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("BLAST public database novelty summary\n")
        fh.write("====================================\n\n")
        fh.write("Hit criteria\n")
        fh.write("------------\n")
        fh.write(f"- evalue <= {args.hit_evalue}\n")
        fh.write(f"- qcovs >= {args.hit_qcovs}\n")
        fh.write(f"- pident >= {args.hit_pident}\n")
        fh.write(
            f"- alignment_length >= max({args.hit_length}, ceil({args.hit_min_aln_frac} * query_length))\n"
        )
        fh.write(f"- db={args.db}, task={args.task}, max_target_seqs={args.max_target_seqs}, max_hsps={args.max_hsps}\n")
        fh.write(f"- remote={'no' if args.no_remote else 'yes'}\n\n")

        for r in rows:
            total = int(r["total_unitigs"])
            unm = int(r["unmapped_unitigs_ref"])
            fh.write(f"{r['dataset']}\n")
            fh.write("-" * len(str(r["dataset"])) + "\n")
            fh.write(f"total_trait_associated_unitigs               : {total}\n")
            fh.write(f"mapped_in_reference                          : {r['mapped_unitigs_ref']} ({pct(int(r['mapped_unitigs_ref']), total)})\n")
            fh.write(f"unmapped_in_reference                        : {unm} ({pct(unm, total)})\n")
            fh.write(f"queries_with_any_blast_row                   : {r['queries_with_any_blast_row']} ({pct(int(r['queries_with_any_blast_row']), total)})\n")
            fh.write(f"queries_with_significant_public_db_hit       : {r['queries_with_significant_hit']} ({pct(int(r['queries_with_significant_hit']), total)})\n")
            fh.write(f"queries_without_significant_hit_novel        : {r['queries_without_significant_hit_novel']} ({pct(int(r['queries_without_significant_hit_novel']), total)})\n")
            fh.write(f"unmapped_with_significant_public_db_hit      : {r['unmapped_with_significant_hit']} ({pct(int(r['unmapped_with_significant_hit']), unm)})\n")
            fh.write(f"unmapped_without_significant_hit_novel       : {r['unmapped_novel']} ({pct(int(r['unmapped_novel']), unm)})\n")
            fh.write(f"blast_rows_total                             : {r['blast_rows_total']}\n")
            fh.write(f"batches                                      : {r['num_batches']} (batch_size={r['batch_size']})\n")
            fh.write(f"batch_progress_log                           : {r['batch_progress_log']}\n")
            fh.write(f"resume_checkpoint_json                       : {r['resume_checkpoint_json']}\n")
            fh.write(f"combined_hits_all_tsv                        : {r['combined_hits_all_tsv']}\n")
            fh.write(f"combined_hits_significant_tsv                : {r['combined_hits_significant_tsv']}\n")
            fh.write(f"result_dir                                   : {r['result_dir']}\n\n")


def write_supplementary5_tsv(path: Path, args: argparse.Namespace, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    def pct_fmt(n: int, d: int) -> str:
        if d <= 0:
            return "n/a"
        return f"{(100.0 * n / d):.2f}%"

    with path.open("w") as fh:
        fh.write("Supplementary File 5\n")
        fh.write("Human-readable BLAST novelty classification summary\n\n")
        fh.write("Counts and percentages are reported for reference mapping status and public DB BLAST support.\n")
        fh.write("This report is dataset-centric and publication-ready (no run-path metadata).\n\n")

        for r in rows:
            dataset = str(r["dataset"])
            total = int(r["total_unitigs"])
            mapped = int(r["mapped_unitigs_ref"])
            unmapped = int(r["unmapped_unitigs_ref"])
            any_hit = int(r["queries_with_any_blast_row"])
            sig_hit = int(r["queries_with_significant_hit"])
            novel = int(r["queries_without_significant_hit_novel"])
            unmapped_hit = int(r["unmapped_with_significant_hit"])
            unmapped_novel = int(r["unmapped_novel"])

            mapped_hit = max(0, sig_hit - unmapped_hit)
            mapped_novel = max(0, mapped - mapped_hit)

            fh.write("==============================================================================\n")
            fh.write(f"{dataset.upper()} -- BLAST NOVELTY SUMMARY\n")
            fh.write("==============================================================================\n\n")

            fh.write("Reference mapping\n\n")
            fh.write(f"{'metric':<58} {'count':>10} {'% total':>12}\n")
            fh.write(f"{'-' * 58} {'-' * 10:>10} {'-' * 12:>12}\n")
            fh.write(f"{'Mapped to available reference':<58} {mapped:>10d} {pct_fmt(mapped, total):>12}\n")
            fh.write(f"{'Not mapped to available reference':<58} {unmapped:>10d} {pct_fmt(unmapped, total):>12}\n\n")

            fh.write("Public database BLAST support\n\n")
            fh.write(f"{'metric':<58} {'count':>10} {'% total':>12}\n")
            fh.write(f"{'-' * 58} {'-' * 10:>10} {'-' * 12:>12}\n")
            fh.write(f"{'Significant public DB BLAST hit':<58} {sig_hit:>10d} {pct_fmt(sig_hit, total):>12}\n")
            fh.write(f"{'No significant public DB BLAST hit':<58} {novel:>10d} {pct_fmt(novel, total):>12}\n\n")

            fh.write("Combined reference and public DB classification\n\n")
            fh.write(f"{'class':<58} {'count':>10} {'% total':>12}\n")
            fh.write(f"{'-' * 58} {'-' * 10:>10} {'-' * 12:>12}\n")
            fh.write(f"{'Mapped reference + significant BLAST hit':<58} {mapped_hit:>10d} {pct_fmt(mapped_hit, total):>12}\n")
            fh.write(f"{'Mapped reference + no significant BLAST hit':<58} {mapped_novel:>10d} {pct_fmt(mapped_novel, total):>12}\n")
            fh.write(f"{'Unmapped reference + significant BLAST hit':<58} {unmapped_hit:>10d} {pct_fmt(unmapped_hit, total):>12}\n")
            fh.write(f"{'Unmapped reference + no significant BLAST hit':<58} {unmapped_novel:>10d} {pct_fmt(unmapped_novel, total):>12}\n\n")

            fh.write("Reference-unmapped subset\n\n")
            fh.write(f"{'class':<58} {'count':>10} {'% unmapped':>14} {'% total':>12}\n")
            fh.write(f"{'-' * 58} {'-' * 10:>10} {'-' * 14:>14} {'-' * 12:>12}\n")
            fh.write(
                f"{'Unmapped reference + significant BLAST hit':<58} "
                f"{unmapped_hit:>10d} {pct_fmt(unmapped_hit, unmapped):>14} {pct_fmt(unmapped_hit, total):>12}\n"
            )
            fh.write(
                f"{'Unmapped reference + no significant BLAST hit':<58} "
                f"{unmapped_novel:>10d} {pct_fmt(unmapped_novel, unmapped):>14} {pct_fmt(unmapped_novel, total):>12}\n"
            )
            fh.write(
                f"{'Total reference-unmapped unitigs':<58} "
                f"{unmapped:>10d} {pct_fmt(unmapped, unmapped):>14} {pct_fmt(unmapped, total):>12}\n\n"
            )

            fh.write("Context across all trait-associated unitigs\n\n")
            fh.write(f"{'metric':<58} {'count':>10} {'% total':>12}\n")
            fh.write(f"{'-' * 58} {'-' * 10:>10} {'-' * 12:>12}\n")
            fh.write(f"{'Total trait-associated unitigs':<58} {total:>10d} {pct_fmt(total, total):>12}\n")
            fh.write(f"{'Mapped to available reference':<58} {mapped:>10d} {pct_fmt(mapped, total):>12}\n")
            fh.write(f"{'Reference-unmapped unitigs':<58} {unmapped:>10d} {pct_fmt(unmapped, total):>12}\n")
            fh.write(f"{'Queries with any BLAST row':<58} {any_hit:>10d} {pct_fmt(any_hit, total):>12}\n")
            fh.write(f"{'Queries with significant BLAST hit (all unitigs)':<58} {sig_hit:>10d} {pct_fmt(sig_hit, total):>12}\n")
            fh.write(f"{'Queries without significant BLAST hit (all unitigs)':<58} {novel:>10d} {pct_fmt(novel, total):>12}\n")
            fh.write(f"{'Total BLAST rows retained':<58} {int(r['blast_rows_total']):>10d} {'n/a':>12}\n\n")

        fh.write("BLAST significance criteria used for all datasets\n\n")
        fh.write(f"{'criterion':<36} {'value':<20}\n")
        fh.write(f"{'-' * 36} {'-' * 20}\n")
        fh.write(f"{'evalue <=':<36} {args.hit_evalue:<20}\n")
        fh.write(f"{'qcovs >=':<36} {args.hit_qcovs:<20}\n")
        fh.write(f"{'pident >=':<36} {args.hit_pident:<20}\n")
        fh.write(f"{'alignment length absolute min (bp)':<36} {args.hit_length:<20}\n")
        fh.write(f"{'alignment length min fraction':<36} {args.hit_min_aln_frac:<20}\n")
        fh.write(f"{'db / task':<36} {args.db} / {args.task}\n")
        fh.write(f"{'max_target_seqs / max_hsps':<36} {args.max_target_seqs} / {args.max_hsps}\n")
        fh.write(f"{'remote mode':<36} {'no' if args.no_remote else 'yes'}\n")


def main() -> int:
    args = parse_args()
    if not args.parse_only and shutil.which(args.blastn_bin) is None:
        print(f"[ERROR] blastn binary not found: {args.blastn_bin}", file=sys.stderr)
        return 2

    datasets_all = [
        DatasetConfig("populus", "Populus tremula", args.populus_fasta, args.populus_unmapped_ids),
        DatasetConfig("ginkgo", "Ginkgo biloba", args.ginkgo_fasta, args.ginkgo_unmapped_ids),
    ]
    selected = set(args.datasets)
    datasets = [ds for ds in datasets_all if ds.key in selected]

    for ds in datasets:
        if not ds.fasta.exists():
            print(f"[ERROR] FASTA not found: {ds.fasta}", file=sys.stderr)
            return 2
        if not ds.unmapped_ids.exists():
            print(f"[ERROR] unmapped ID file not found: {ds.unmapped_ids}", file=sys.stderr)
            return 2

    results = []
    for ds in datasets:
        log(f"Processing {ds.label}")
        results.append(process_dataset(args, ds))

    write_summary_txt(args.summary_txt, args, results)
    write_supplementary5_tsv(args.supplementary5_tsv, args, results)
    log(f"Wrote summary: {args.summary_txt}")
    log(f"Wrote Supplementary File 5 TSV: {args.supplementary5_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
