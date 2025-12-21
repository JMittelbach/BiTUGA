import sys
sys.dont_write_bytecode = True

import argparse
from pathlib import Path
from datetime import datetime
import math


class Logger:
    def __init__(self, logging_dir: Path, pass_name: str = "unitigs_fisher"):
        try:
            logging_dir.mkdir(parents=True, exist_ok=True)
            self.main_log = (logging_dir / "main.log").open("a")
            self.pass_log = (logging_dir / f"{pass_name}.log").open("a")
            self.use_files = True
        except Exception:
            self.main_log = sys.stdout
            self.pass_log = sys.stdout
            self.use_files = False

    def log(self, msg: str, level: str = "INFO"):
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{ts}] [{level}] {msg}\n"
        self.main_log.write(line)
        self.pass_log.write(line)
        self.main_log.flush()
        self.pass_log.flush()

    def close(self):
        if self.use_files:
            self.main_log.close()
            if self.pass_log is not self.main_log:
                self.pass_log.close()


def parse_stats(detailed_path: Path, min_output_delta: float, use_fdr: bool, alpha: float, p_threshold: float):
    with detailed_path.open("r") as fh:
        header = fh.readline().rstrip("\n")
        if not header:
            raise SystemExit("unitig_association_details.tsv is empty")
        cols = header.split("\t")
        col_idx = {name: i for i, name in enumerate(cols)}

        prev_cols = [name for name in cols if name.startswith("prev_")]
        if len(prev_cols) != 2:
            raise SystemExit(f"expected exactly two prev_* columns, found: {prev_cols}")
        prev1_col = prev_cols[0]
        prev2_col = prev_cols[1]
        trait1 = prev1_col[len("prev_"):]
        trait2 = prev2_col[len("prev_"):]

        count1_col = "count_" + trait1
        count2_col = "count_" + trait2
        if count1_col not in col_idx or count2_col not in col_idx:
            raise SystemExit(f"expected count columns {count1_col} and {count2_col} in detailed stats")
        if "delta_prev" not in col_idx:
            raise SystemExit("delta_prev column not found in detailed stats")

        idx_unitig = col_idx["unitig"]
        idx_delta = col_idx["delta_prev"]
        idx_count1 = col_idx[count1_col]
        idx_count2 = col_idx[count2_col]
        idx_bias = col_idx["bias"]
        use_storey = "q_value" in col_idx
        fdr_available = "q_value" in col_idx or "p_adj" in col_idx
        if use_fdr and fdr_available:
            adj_col = "q_value" if use_storey else "p_adj"
            idx_adj = col_idx[adj_col]
            use_fdr_final = True
        else:
            use_fdr_final = False
            if "p_value" not in col_idx:
                raise SystemExit("p_value column not found in detailed stats")
            idx_p = col_idx["p_value"]
            p_threshold = p_threshold if p_threshold is not None else 0.01

        significant = {}
        by_bias_t1 = set()
        by_bias_t2 = set()
        stat_label = "q" if (use_storey and use_fdr_final) else ("p_adj" if use_fdr_final else "p")

        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            uid = parts[idx_unitig]
            bias = parts[idx_bias]
            if bias == trait1:
                bias_key = "trait1"
            elif bias == trait2:
                bias_key = "trait2"
            elif bias == "trait1" or bias == "trait2":
                bias_key = bias
            else:
                bias_key = bias
            if use_fdr_final:
                p_adj_str = parts[idx_adj]
                if p_adj_str == "NA":
                    continue
                try:
                    q = float(p_adj_str)
                except ValueError:
                    continue
            else:
                p_raw_str = parts[idx_p]
                if p_raw_str == "NA":
                    continue
                try:
                    q = float(p_raw_str)
                except ValueError:
                    continue

            if use_fdr_final:
                is_sig = q <= alpha
            else:
                is_sig = q <= p_threshold
            if not is_sig:
                continue

            delta_str = parts[idx_delta]
            if delta_str == "NA":
                continue
            try:
                delta = float(delta_str)
            except ValueError:
                continue
            if delta < min_output_delta:
                continue

            count1 = int(parts[idx_count1])
            count2 = int(parts[idx_count2])

            significant[uid] = (count1, count2, q, bias)
            if bias_key == "trait1":
                by_bias_t1.add(uid)
            elif bias_key == "trait2":
                by_bias_t2.add(uid)

        return significant, by_bias_t1, by_bias_t2, trait1, trait2, stat_label


def iter_fasta(path: Path):
    with path.open("r") as fh:
        header = None
        seq_lines = []
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line.strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            yield header, "".join(seq_lines)


def main():
    parser = argparse.ArgumentParser(description="PASS VI: extract significant unitigs as FASTA.")
    parser.add_argument("--outdir", default=None, help="Output directory root; default is the parent of this script (project root).")
    parser.add_argument("--unitig-fasta", default=None, help="Path to unitigs.fa (default: <root>/results/unitigs.fa).")
    parser.add_argument("--stats-detailed", default=None, help="Path to unitig_association_details.tsv (default: <root>/results/statistics/unitig_association_details.tsv; legacy: unitig_fishers_detailed.tsv).")
    parser.add_argument("--min-bias-prev", type=float, default=0.0, help="Minimum required delta_prev for FASTA output (default: 0.0).")
    parser.add_argument("--alpha", type=float, default=0.05, help="FDR level when BH is used (default: 0.05).")
    parser.add_argument("--no-fdr", action="store_true", help="Disable FDR correction and filter by raw p-values (stored in p_adj).")
    parser.add_argument("--p-threshold", type=float, default=0.05, help="p-value threshold used when --no-fdr is set (default: 0.05).")
    parser.add_argument("--split-by-bias", action="store_true", help="Write separate FASTA files per biased trait.")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    if args.outdir:
        root = Path(args.outdir).resolve()
    else:
        root = project_root

    results_dir = root / "results"
    stats_dir = results_dir / "statistics"
    logging_dir = root / "logging"
    results_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)
    logger = Logger(logging_dir)

    use_fdr = not args.no_fdr

    if args.unitig_fasta:
        fasta_path = Path(args.unitig_fasta)
    else:
        fasta_path = results_dir / "unitigs.fa"
    if not fasta_path.is_file():
        logger.log(f"unitig fasta not found: {fasta_path}", level="ERROR")
        logger.close()
        raise SystemExit(1)

    if args.stats_detailed:
        detailed_path = Path(args.stats_detailed)
    else:
        detailed_path = stats_dir / "unitig_association_details.tsv"
    if not detailed_path.is_file():
        logger.log(f"detailed stats not found: {detailed_path}", level="ERROR")
        logger.close()
        raise SystemExit(1)
    if not detailed_path.is_file():
        legacy = stats_dir / "unitigs_fishers_detailed.tsv"
        if legacy.is_file():
            detailed_path = legacy
        else:
            logger.log(f"detailed stats not found: {detailed_path}", level="ERROR")
            logger.close()
            raise SystemExit(1)

    logger.log(
        "Extract arguments: "
        f"root={root.name}, fasta={fasta_path.name}, stats_detailed={detailed_path.name}, "
        f"min_output_delta={args.min_bias_prev}, "
        f"use_fdr={use_fdr}, alpha={args.alpha}, p_threshold={args.p_threshold}, "
        f"split_by_bias={args.split_by_bias}"
    )

    p_threshold_val = args.p_threshold if args.p_threshold is not None else 0.01

    significant, by_bias_t1, by_bias_t2, trait1, trait2, stat_label = parse_stats(
        detailed_path=detailed_path,
        min_output_delta=args.min_bias_prev,
        use_fdr=use_fdr,
        alpha=args.alpha,
        p_threshold=p_threshold_val,
    )

    n_sig = len(significant)
    logger.log(f"Significant unitigs after output delta_prev filter: {n_sig}")
    logger.log(f"Biased towards {trait1}: {len(by_bias_t1)}")
    logger.log(f"Biased towards {trait2}: {len(by_bias_t2)}")

    if n_sig == 0:
        out_global = results_dir / "significant_unitigs.fa"
        out_global.open("w").close()
        logger.log(f"No significant unitigs; wrote empty FASTA to {out_global}")
        logger.close()
        return

    out_global = results_dir / "significant_unitigs.fa"
    out_t1 = results_dir / f"significant_unitigs_{trait1}.fa"
    out_t2 = results_dir / f"significant_unitigs_{trait2}.fa"

    sorted_uids = sorted(significant.items(), key=lambda kv: kv[1][2])

    seq_map = {}
    needed = set(significant.keys())
    for header, seq in iter_fasta(fasta_path):
        uid = header[1:].split()[0]
        if uid in needed:
            seq_map[uid] = seq

    f_global = out_global.open("w")
    f_t1 = out_t1.open("w") if args.split_by_bias else None
    f_t2 = out_t2.open("w") if args.split_by_bias else None

    written_global = 0
    written_t1 = 0
    written_t2 = 0

    for uid, (count1, count2, q, bias) in sorted_uids:
        seq = seq_map.get(uid)
        if seq is None:
            continue
        q_str = f"{q:.6g}"
        c1_str = str(count1)
        c2_str = str(count2)
        new_header = f">{uid} {trait1}={c1_str} {trait2}={c2_str} {stat_label}={q_str}"
        f_global.write(new_header + "\n")
        f_global.write(seq + "\n")
        written_global += 1

        if args.split_by_bias:
            if bias == "trait1" and f_t1 is not None:
                f_t1.write(new_header + "\n")
                f_t1.write(seq + "\n")
                written_t1 += 1
            elif bias == "trait2" and f_t2 is not None:
                f_t2.write(new_header + "\n")
                f_t2.write(seq + "\n")
                written_t2 += 1

    f_global.close()
    if f_t1 is not None:
        f_t1.close()
    if f_t2 is not None:
        f_t2.close()

    top_n = min(10, len(sorted_uids))
    if top_n > 0:
        logger.log(f"Top {top_n} significant unitigs ({stat_label} ascending):")
        for uid, (count1, count2, q, bias) in sorted_uids[:top_n]:
            seq = seq_map.get(uid, "")
            logger.log(f"  {uid}\t{stat_label}={q:.6g}\t{trait1}={count1}\t{trait2}={count2}\tbias={bias}\tseq={seq}")

        logger.log(f"Wrote {written_global} unitigs to {out_global.name}")
        if args.split_by_bias:
            logger.log(f"Wrote {written_t1} unitigs biased towards {trait1} to {out_t1.name}")
            logger.log(f"Wrote {written_t2} unitigs biased towards {trait2} to {out_t2.name}")
    logger.close()


if __name__ == "__main__":
    main()
