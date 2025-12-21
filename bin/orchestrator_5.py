import sys
sys.dont_write_bytecode = True

import argparse
import subprocess
import gzip
from pathlib import Path
from datetime import datetime
from collections import Counter
from typing import Dict, Iterable, List


def log_stdout(msg: str):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    sys.stdout.write(f"[{ts}] [INFO] {msg}\n")
    sys.stdout.flush()


def run_cmd(cmd, log_cmd=None, cwd=None):
    pretty = " ".join(str(x) for x in (log_cmd if log_cmd is not None else cmd))
    log_stdout("Running: " + pretty)
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        raise SystemExit(result.returncode)

def load_unitig_lengths(fa_path: Path) -> Dict[str, int]:
    lens: Dict[str, int] = {}
    if not fa_path.is_file():
        return lens
    opener = gzip.open if str(fa_path).endswith(".gz") else open
    with opener(fa_path, "rt") as fh:
        uid = None
        seq_parts: List[str] = []
        for line in fh:
            if line.startswith(">"):
                if uid is not None:
                    lens[uid] = len("".join(seq_parts))
                uid = line[1:].strip().split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if uid is not None:
            lens[uid] = len("".join(seq_parts))
    return lens

def lengths_from_fasta(fa_path: Path, len_map: Dict[str, int]) -> List[int]:
    lengths: List[int] = []
    if not fa_path.is_file():
        return lengths
    opener = gzip.open if str(fa_path).endswith(".gz") else open
    with opener(fa_path, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            uid = line[1:].strip().split()[0]
            if uid in len_map:
                lengths.append(len_map[uid])
    return lengths

def append_histogram(lengths: Iterable[int], dest: Path, label: str):
    counts = Counter(lengths)
    if not counts:
        return
    with dest.open("a", encoding="utf-8") as fh:
        fh.write(f"# {label} length histogram (len\tcount)\n")
        for L in sorted(counts):
            fh.write(f"{L}\t{counts[L]}\n")
        fh.write("\n")


def main():
    parser = argparse.ArgumentParser(description="PASS V orchestrator: Fisher tests and unitig extraction.")
    parser.add_argument("--outdir", default=None, help="Output directory root; default is the parent of this script (project root).")
    parser.add_argument("--metadata", default=None, help="Path to metadata.tsv (default: <project_root>/metadata.tsv).")
    parser.add_argument("--min-global-prev", type=float, default=None, help="Minimum global prevalence to test a unitig (default: 0.10).")
    parser.add_argument("--max-global-prev", type=float, default=None, help="Maximum global prevalence to test a unitig (default: 0.90).")
    parser.add_argument("--adaptive-prev-thresholds-unitig", choices=["auto", "off"], default="auto",
                        help="Adaptive prevalence defaults for asymmetric groups (auto) or off for static.")
    parser.add_argument("--min-delta-prev", type=float, default=0.0, help="Minimum absolute prevalence difference to test a unitig (default: 0.0).")
    parser.add_argument("--min-bias-prev", type=float, default=0.0, help="Minimum prevalence in biased trait for FASTA output (default: 0.0).")
    parser.add_argument("--alpha", type=float, default=0.05, help="FDR level for BH correction (default: 0.05).")
    parser.add_argument("--sig-mode", dest="fdr_mode", help=argparse.SUPPRESS)
    parser.add_argument("--fdr", dest="fdr_mode", choices=["auto", "bh", "storey", "none"], default="auto",
                        help="Multiple testing: auto (default) chooses BH if feasible, otherwise raw p; 'bh' forces BH; 'storey' uses Storey q-values; 'none' uses raw p threshold.")
    parser.add_argument("--p-threshold", type=float, default=0.01, help="p-value threshold used when fdr-mode=none or auto falls back (default: 0.01).")
    parser.add_argument("--hist-bins", type=int, default=50, help="Number of bins for p-value and p_adj histograms (default: 50).")
    parser.add_argument("--split-by-bias", action="store_true", help="Write separate FASTA files per biased trait.")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    if args.outdir:
        root = Path(args.outdir).resolve()
    else:
        root = project_root

    results_dir = root / "results"
    logging_dir = root / "logging"
    results_dir.mkdir(parents=True, exist_ok=True)
    logging_dir.mkdir(parents=True, exist_ok=True)

    log_stdout(f"PASS V orchestrator starting with root={root.name}")

    fisher_script = script_dir / "unitig_fisher_5.py"
    extract_script = script_dir / "extract_unitigs_5.py"

    if args.metadata:
        metadata_path = Path(args.metadata)
    else:
        metadata_path = project_root / "metadata.tsv"

    fisher_cmd_log = [
        "python3", fisher_script.name,
        "--outdir", ".",
        "--metadata", metadata_path.name,
        "--adaptive-prev-thresholds-unitig", str(args.adaptive_prev_thresholds_unitig),
        "--min-delta-prev", str(args.min_delta_prev),
        "--alpha", str(args.alpha),
        "--fdr", str(args.fdr_mode),
        "--hist-bins", str(args.hist_bins),
    ]
    fisher_cmd_exec = [
        sys.executable, str(fisher_script),
        "--outdir", str(root),
        "--metadata", str(metadata_path),
        "--adaptive-prev-thresholds-unitig", str(args.adaptive_prev_thresholds_unitig),
        "--min-delta-prev", str(args.min_delta_prev),
        "--alpha", str(args.alpha),
        "--fdr", str(args.fdr_mode),
        "--hist-bins", str(args.hist_bins),
    ]
    if args.min_global_prev is not None:
        fisher_cmd_log.extend(["--min-global-prev", str(args.min_global_prev)])
        fisher_cmd_exec.extend(["--min-global-prev", str(args.min_global_prev)])
    if args.max_global_prev is not None:
        fisher_cmd_log.extend(["--max-global-prev", str(args.max_global_prev)])
        fisher_cmd_exec.extend(["--max-global-prev", str(args.max_global_prev)])
    fisher_cmd_log.extend(["--p-threshold", f"{args.p_threshold}"])
    fisher_cmd_exec.extend(["--p-threshold", f"{args.p_threshold}"])

    run_cmd(fisher_cmd_exec, log_cmd=fisher_cmd_log, cwd=root)

    extract_cmd_log = [
        "python3", extract_script.name,
        "--outdir", ".",
        "--min-bias-prev", str(args.min_bias_prev),
        "--alpha", str(args.alpha),
        "--p-threshold", f"{args.p_threshold}",
    ]
    extract_cmd_exec = [
        sys.executable, str(extract_script),
        "--outdir", str(root),
        "--min-bias-prev", str(args.min_bias_prev),
        "--alpha", str(args.alpha),
        "--p-threshold", f"{args.p_threshold}",
    ]
    if args.fdr_mode == "none":
        extract_cmd_log.append("--no-fdr")
        extract_cmd_exec.append("--no-fdr")
    if args.split_by_bias:
        extract_cmd_log.append("--split-by-bias")
        extract_cmd_exec.append("--split-by-bias")

    run_cmd(extract_cmd_exec, log_cmd=extract_cmd_log, cwd=root)

    try:
        unitigs_fa = results_dir / "unitigs.fa"
        sig_fa = results_dir / "significant_unitigs.fa"
        fishers_txt = results_dir / "statistics" / "unitig_association_details.tsv"
        len_map = load_unitig_lengths(unitigs_fa)
        sig_lengths = lengths_from_fasta(sig_fa, len_map)
        append_histogram(sig_lengths, fishers_txt, "significant unitigs")
    except Exception as exc:
        log_stdout(f"[WARN] Failed to append significant unitig length histogram: {exc}")

    log_stdout("PASS V orchestrator finished")


if __name__ == "__main__":
    main()
