import sys
sys.dont_write_bytecode = True

import argparse
import math
from pathlib import Path
from typing import Dict, Tuple, List
from datetime import datetime
from fdr_5 import benjamini_hochberg, storey_qvalues


class Logger:
    def __init__(self, logging_dir: Path, pass_name: str = "unitigs_fisher"):
        logging_dir.mkdir(parents=True, exist_ok=True)
        self.main_log = (logging_dir / "main.log").open("a")
        self.pass_log = (logging_dir / f"{pass_name}.log").open("a")

    def log(self, msg: str, level: str = "INFO"):
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{ts}] [{level}] {msg}\n"
        self.main_log.write(line)
        self.pass_log.write(line)
        self.main_log.flush()
        self.pass_log.flush()

    def close(self):
        self.main_log.close()
        self.pass_log.close()


def fisher_two_sided(a: int, b: int, c: int, d: int) -> float:
    from math import comb
    r1 = a + b
    r2 = c + d
    col1 = a + c
    N = r1 + r2
    x_min = max(0, col1 - r2)
    x_max = min(col1, r1)
    denom = comb(N, col1)

    def prob(x: int) -> float:
        return comb(r1, x) * comb(r2, col1 - x) / denom

    p_obs = prob(a)
    p_sum = 0.0
    eps = 1e-12
    for x in range(x_min, x_max + 1):
        p = prob(x)
        if p <= p_obs + eps:
            p_sum += p
    if p_sum > 1.0:
        p_sum = 1.0
    return p_sum


def build_fisher_lookup(n1: int, n2: int) -> Dict[Tuple[int, int], float]:
    lookup: Dict[Tuple[int, int], float] = {}
    for a in range(n1 + 1):
        for c in range(n2 + 1):
            b = n1 - a
            d = n2 - c
            p = fisher_two_sided(a, b, c, d)
            lookup[(a, c)] = p
    return lookup


def read_metadata(meta_path: Path, trait1_arg: str, trait2_arg: str):
    with meta_path.open("r") as fh:
        header = fh.readline().rstrip("\n")
        if not header:
            raise SystemExit("metadata file is empty")
        cols = header.split("\t")
        idx_sample = cols.index("sample_id")
        idx_trait = cols.index("trait")
        sample_to_trait = {}
        traits = set()
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            sid = parts[idx_sample]
            tr = parts[idx_trait]
            if sid in sample_to_trait and sample_to_trait[sid] != tr:
                raise SystemExit(f"sample_id {sid} has inconsistent traits in metadata")
            sample_to_trait[sid] = tr
            traits.add(tr)
    if trait1_arg and trait2_arg:
        trait1 = trait1_arg
        trait2 = trait2_arg
    else:
        if len(traits) != 2:
            raise SystemExit(f"expected exactly 2 distinct traits in metadata, found {len(traits)}: {sorted(traits)}")
        trait1, trait2 = sorted(traits)
    n1 = sum(1 for t in sample_to_trait.values() if t == trait1)
    n2 = sum(1 for t in sample_to_trait.values() if t == trait2)
    total = len(sample_to_trait)
    if n1 == 0 or n2 == 0:
        raise SystemExit(f"one of the traits has zero individuals: {trait1}={n1}, {trait2}={n2}")
    return trait1, trait2, n1, n2, total


def main():
    parser = argparse.ArgumentParser(description="PASS V: Fisher tests on unitig prevalence counts.")
    parser.add_argument("--outdir", default=None, help="Output directory root; default is the parent of this script (project root).")
    parser.add_argument("--input", help="Optional explicit path to unitig_presence.tsv (overrides default).")
    parser.add_argument("--metadata", default=None, help="Path to metadata.tsv (default: <project_root>/metadata.tsv).")
    parser.add_argument("--trait1", default=None, help="Name of first trait; if not set, inferred from metadata.")
    parser.add_argument("--trait2", default=None, help="Name of second trait; if not set, inferred from metadata.")
    parser.add_argument("--min-global-prev", type=float, default=None, help="Minimum global prevalence to test a unitig (default: 0.10).")
    parser.add_argument("--max-global-prev", type=float, default=None, help="Maximum global prevalence to test a unitig (default: 0.90).")
    parser.add_argument("--adaptive-prev-thresholds-unitig", dest="adaptive_prev_thresholds", choices=["auto", "off"], default="auto",
                        help="auto adapts defaults for asymmetric groups; off keeps static thresholds (defaults 0.10/0.90).")
    parser.add_argument("--adaptive-prev-thresholds", dest="adaptive_prev_thresholds", choices=["auto", "off"], help=argparse.SUPPRESS)
    parser.add_argument("--min-delta-prev", type=float, default=0.0, help="Minimum absolute prevalence difference to test a unitig (default: 0.0).")
    parser.add_argument("--alpha", type=float, default=0.05, help="FDR level for BH correction (default: 0.05).")
    parser.add_argument("--fdr", choices=["auto", "bh", "storey", "none"], default="auto",
                        help="Multiple testing: auto (default) chooses BH if feasible, otherwise raw p; 'bh' forces BH; 'storey' uses Storey q-values; 'none' uses raw p threshold.")
    parser.add_argument("--sig-mode", dest="fdr", help=argparse.SUPPRESS)
    parser.add_argument("--p-threshold", type=float, default=0.01, help="p-value threshold used when fdr=none or auto falls back (default: 0.01).")
    parser.add_argument("--hist-bins", type=int, default=200, help="Number of bins for p-value and adjusted/q-value histograms (default: 200).")
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

    mode = args.fdr

    if args.metadata:
        meta_path = Path(args.metadata)
    else:
        meta_path = project_root / "metadata.tsv"
    if not meta_path.is_file():
        logger.log(f"metadata file not found: {meta_path}", level="ERROR")
        logger.close()
        raise SystemExit(1)

    trait1, trait2, n1, n2, total_n = read_metadata(meta_path, args.trait1, args.trait2)

    min_set = args.min_global_prev is not None
    max_set = args.max_global_prev is not None
    base_min = args.min_global_prev if min_set else 0.10
    base_max = args.max_global_prev if max_set else 0.90

    minority_ratio = (min(n1, n2) / total_n) if total_n > 0 else 0.0
    adaptive_min = minority_ratio * 0.5
    adaptive_max = 1.0 - adaptive_min

    adaptive_applied = args.adaptive_prev_thresholds == "auto" and not min_set and not max_set and n1 != n2
    eff_min = adaptive_min if adaptive_applied else base_min
    eff_max = adaptive_max if adaptive_applied else base_max

    min_carriers = max(2, math.ceil(eff_min * total_n))
    max_carriers = max(0, math.floor(eff_max * total_n))

    def rel(p: Path) -> str:
        try:
            return "." + str(p.resolve().absolute()).replace(str(root), "")
        except Exception:
            return str(p)

    logger.log(
        "Arguments: "
        f"project_root={project_root.name}, root={root.name}, outdir_arg={args.outdir if not args.outdir else Path(args.outdir).name}, "
        f"input={args.input if not args.input else Path(args.input).name}, metadata={meta_path.name}, "
        f"trait1_arg={args.trait1}, trait2_arg={args.trait2}, "
        f"min_global_prev={eff_min}, "
        f"max_global_prev={eff_max}, "
        f"adaptive_prev_mode={args.adaptive_prev_thresholds}, adaptive_prev_applied={adaptive_applied}, "
        f"minority_ratio={minority_ratio:.4f}, "
        f"min_global_prev_min_carriers={min_carriers}, "
        f"max_global_prev_max_carriers={max_carriers}, "
        f"min_delta_prev={args.min_delta_prev}, "
        f"alpha={args.alpha}, fdr_mode={mode}, p_threshold={args.p_threshold}, "
        f"hist_bins={args.hist_bins}"
    )
    logger.log(f"Traits: trait1={trait1} (n={n1}), trait2={trait2} (n={n2}), total={total_n}")

    if args.input:
        presence_path = Path(args.input)
    else:
        presence_path = results_dir / "unitig_presence.tsv"
    if not presence_path.is_file():
        logger.log(f"input file not found: {presence_path}", level="ERROR")
        logger.close()
        raise SystemExit(1)

    logger.log("Building Fisher lookup table")
    lookup = build_fisher_lookup(n1, n2)
    logger.log(f"Lookup size: {len(lookup)} combinations")

    stats_rows: List[dict] = []
    tested_indices: List[int] = []
    p_values: List[float] = []

    filtered_by_prev_threshold_min = 0
    filtered_by_prev_threshold_max = 0
    filtered_by_prev_threshold_both = 0
    filtered_min_delta_prev = 0

    pre_bias_t1 = 0
    pre_bias_t2 = 0
    pre_bias_tie = 0

    with presence_path.open("r") as fin:
        header = fin.readline()
        if not header:
            logger.log("unitig_presence.tsv is empty", level="ERROR")
            logger.close()
            raise SystemExit(1)
        line_no = 1
        for line in fin:
            line_no += 1
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:
                logger.log(f"skipping malformed line {line_no} in {presence_path}", level="WARN")
                continue
            unitig_id = parts[0]
            try:
                a = int(parts[1])
                c = int(parts[2])
            except ValueError:
                logger.log(f"non-integer counts in line {line_no}, skipping", level="WARN")
                continue

            total_carriers = a + c
            prev1 = a / n1 if n1 > 0 else math.nan
            prev2 = c / n2 if n2 > 0 else math.nan
            global_prev = total_carriers / total_n if total_n > 0 else math.nan
            delta_prev = abs(prev1 - prev2) if not math.isnan(prev1) and not math.isnan(prev2) else math.nan

            if math.isnan(prev1) or math.isnan(prev2):
                bias = "tie"
                pre_bias_tie += 1
            elif prev1 > prev2:
                bias = "trait1"
                pre_bias_t1 += 1
            elif prev2 > prev1:
                bias = "trait2"
                pre_bias_t2 += 1
            else:
                bias = "tie"
                pre_bias_tie += 1

            b = n1 - a
            d = n2 - c
            a_c = a + 0.5
            b_c = b + 0.5
            c_c = c + 0.5
            d_c = d + 0.5
            odds_ratio = (a_c * d_c) / (b_c * c_c)

            fails_prev = total_carriers < min_carriers
            fails_prev_max = total_carriers > max_carriers if max_carriers > 0 else False
            fails_min_delta = delta_prev < args.min_delta_prev

            if fails_prev:
                filtered_by_prev_threshold_min += 1
            if fails_prev_max:
                filtered_by_prev_threshold_max += 1
            if fails_prev and fails_prev_max:
                filtered_by_prev_threshold_both += 1
            if fails_min_delta:
                filtered_min_delta_prev += 1

            test_flag = not (fails_prev or fails_prev_max or fails_min_delta)

            row = {
                "unitig": unitig_id,
                "count_trait1": a,
                "count_trait2": c,
                "prev_trait1": prev1,
                "prev_trait2": prev2,
                "global_prev": global_prev,
                "delta_prev": delta_prev,
                "odds_ratio": odds_ratio,
                "bias": bias,
                "p_value": math.nan,
                "p_adj": math.nan,
            }
            if test_flag:
                p = lookup[(a, c)]
                row["p_value"] = p
                tested_indices.append(len(stats_rows))
                p_values.append(p)
            stats_rows.append(row)

    n_total = len(stats_rows)
    n_tested = len(tested_indices)
    logger.log(f"Unitigs read: {n_total}")
    logger.log(f"Unitigs tested after filters: {n_tested}")
    logger.log(f"Filtered by prevalence threshold min (min_carriers={min_carriers}): {filtered_by_prev_threshold_min}")
    logger.log(f"Filtered by prevalence threshold max (max_carriers={max_carriers}): {filtered_by_prev_threshold_max}")
    logger.log(f"Filtered by prevalence threshold both (min&&max): {filtered_by_prev_threshold_both}")
    logger.log(f"Filtered by min_delta_prev: {filtered_min_delta_prev}")
    logger.log(f"Pre-bias counts (before filters): trait1={pre_bias_t1}, trait2={pre_bias_t2}, tie={pre_bias_tie}")

    alpha = args.alpha
    p_threshold = args.p_threshold
    p_min_extreme = fisher_two_sided(n1, 0, 0, n2)
    alpha_over_m = alpha / n_tested if n_tested > 0 else math.nan
    bh_feasible = n_tested > 0 and p_min_extreme <= alpha_over_m

    if mode == "bh":
        use_fdr = True
        use_storey = False
    elif mode == "storey":
        use_fdr = True
        use_storey = True
    elif mode in ("p", "none"):
        use_fdr = False
        use_storey = False
    else:
        use_fdr = bh_feasible
        use_storey = False

    mode_resolved = "storey" if use_storey else ("BH" if use_fdr else "raw")
    logger.log(f"Multiple testing mode resolved: mode_arg={mode} -> {mode_resolved}; bh_feasible={bh_feasible} (p_min_extreme={p_min_extreme:.3g}, alpha_over_m={alpha_over_m if not math.isnan(alpha_over_m) else 'NA'})")

    adj_label = "q_value" if use_storey else ("p_adj" if use_fdr else None)

    if use_storey:
        for row in stats_rows:
            row["q_value"] = math.nan

    pi0_est = math.nan

    if use_fdr:
        if p_values:
            if use_storey:
                p_adj_list, pi0 = storey_qvalues(p_values)
                pi0_est = pi0
                logger.log(f"Storey pi0 estimate: {pi0:.4g}")
            else:
                p_adj_list = benjamini_hochberg(p_values)
            for idx_in_list, row_index in enumerate(tested_indices):
                stats_rows[row_index][adj_label] = p_adj_list[idx_in_list]

    sig_total = 0
    sig_t1 = 0
    sig_t2 = 0
    sig_tie = 0

    for row in stats_rows:
        q = row.get(adj_label, math.nan) if use_fdr else row["p_value"]
        if math.isnan(q):
            continue
        if use_fdr:
            is_sig = q <= alpha
        else:
            is_sig = q <= p_threshold
        if is_sig:
            sig_total += 1
            if row["bias"] == "trait1":
                sig_t1 += 1
            elif row["bias"] == "trait2":
                sig_t2 += 1
            else:
                sig_tie += 1

    mode_label = "Storey" if use_storey else ("BH" if use_fdr else "no_fdr")
    logger.log(f"Significant unitigs (mode={mode_label}): {sig_total}")
    logger.log(f"{trait1}-biased: {sig_t1}")
    logger.log(f"{trait2}-biased: {sig_t2}")
    logger.log(f"tie: {sig_tie}")

    target_bins = max(args.hist_bins, 10)
    if target_bins <= 100:
        fine_bins = max(1, target_bins // 2)
        coarse_bins = max(1, target_bins - fine_bins)
    else:
        fine_bins = 100
        coarse_bins = max(1, target_bins - fine_bins)
    fine_width = 0.1 / fine_bins
    coarse_width = 0.9 / coarse_bins
    bin_count = fine_bins + coarse_bins
    hist_counts_p = [0] * bin_count
    hist_counts_q = [0] * bin_count if use_fdr else None

    def bin_index(val: float) -> int:
        if val < 0.0:
            val = 0.0
        if val >= 1.0:
            return bin_count - 1
        if val <= 0.1:
            idx = int(val / fine_width)
            if idx >= fine_bins:
                idx = fine_bins - 1
            return idx
        idx = int((val - 0.1) / coarse_width)
        if idx >= coarse_bins:
            idx = coarse_bins - 1
        return fine_bins + idx

    for row_idx in tested_indices:
        p = stats_rows[row_idx]["p_value"]
        q = stats_rows[row_idx].get(adj_label, math.nan)

        if not math.isnan(p):
            hist_counts_p[bin_index(p)] += 1

        if use_fdr and not math.isnan(q):
            hist_counts_q[bin_index(q)] += 1

    def sort_key(row):
        if use_fdr:
            q = row.get(adj_label, math.nan)
            return (math.isnan(q), q if not math.isnan(q) else 0.0)
        p = row.get("p_value", math.nan)
        return (math.isnan(p), p if not math.isnan(p) else 0.0)

    stats_rows_sorted = sorted(stats_rows, key=sort_key)

    summary_path = stats_dir / "unitig_association_stats.txt"

    count1_name = f"count_{trait1}"
    count2_name = f"count_{trait2}"
    prev1_name = f"prev_{trait1}"
    prev2_name = f"prev_{trait2}"

    detailed_path = stats_dir / "unitig_association_details.tsv"

    p_min_extreme = fisher_two_sided(n1, 0, 0, n2)
    alpha_over_m = alpha / n_tested if n_tested > 0 else math.nan
    bh_can_be_significant = use_fdr and n_tested > 0 and p_min_extreme <= alpha_over_m

    with summary_path.open("w") as fout:
        fout.write(f"total_unitigs\t{n_total}\n")
        fout.write(f"tested_unitigs\t{n_tested}\n")
        fout.write(f"filtered_by_prev_threshold_min\t{filtered_by_prev_threshold_min}\n")
        fout.write(f"filtered_by_prev_threshold_max\t{filtered_by_prev_threshold_max}\n")
        fout.write(f"filtered_by_min_delta_prev\t{filtered_min_delta_prev}\n")
        fout.write(f"filtered_by_prev_threshold_both\t{filtered_by_prev_threshold_both}\n")
        fout.write(f"pre_bias_{trait1}\t{pre_bias_t1}\n")
        fout.write(f"pre_bias_{trait2}\t{pre_bias_t2}\n")
        fout.write(f"pre_bias_tie\t{pre_bias_tie}\n")
        if use_fdr:
            fout.write(f"multiple_testing_correction\t{'storey' if use_storey else 'BH'}\n")
            fout.write(f"alpha\t{alpha:.6f}\n")
            if use_storey:
                fout.write(f"storey_pi0\t{pi0_est:.6g}\n")
        else:
            fout.write(f"multiple_testing_correction\tnone\n")
            fout.write(f"p_threshold\t{p_threshold:.6f}\n")
        fout.write(f"significant_unitigs\t{sig_total}\n")
        fout.write(f"significant_bias_{trait1}\t{sig_t1}\n")
        fout.write(f"significant_bias_{trait2}\t{sig_t2}\n")
        fout.write(f"significant_bias_tie\t{sig_tie}\n")
        fout.write(f"min_global_prev_fraction\t{eff_min}\n")
        fout.write(f"min_global_prev_min_carriers\t{min_carriers}\n")
        fout.write(f"max_global_prev_fraction\t{eff_max}\n")
        fout.write(f"max_global_prev_max_carriers\t{max_carriers}\n")
        fout.write(f"adaptive_prev_mode\t{args.adaptive_prev_thresholds}\n")
        fout.write(f"adaptive_prev_applied\t{str(adaptive_applied).lower()}\n")
        if adaptive_applied:
            fout.write(f"adaptive_prev_min_fraction\t{adaptive_min}\n")
            fout.write(f"adaptive_prev_max_fraction\t{adaptive_max}\n")
            fout.write(f"minority_ratio\t{minority_ratio}\n")
        fout.write(f"min_delta_prev\t{args.min_delta_prev}\n")
        fout.write(f"metadata\t{meta_path}\n")
        fout.write(f"trait1\t{trait1}\n")
        fout.write(f"trait2\t{trait2}\n")
        fout.write(f"hist_bins\t{bin_count}\n")
        fout.write(f"p_min_extreme\t{p_min_extreme:.6g}\n")
        if not math.isnan(alpha_over_m):
            fout.write(f"alpha_over_m\t{alpha_over_m:.6g}\n")
        else:
            fout.write("alpha_over_m\tNA\n")
        fout.write(f"bh_can_be_significant\t{str(bh_can_be_significant).lower()}\n")
        hist_header = "pvalue_qvalue_histogram" if use_storey else ("pvalue_padj_histogram" if use_fdr else "pvalue_histogram")
        fout.write(f"{hist_header}\n")
        if use_fdr:
            adj_colname = "count_q_value" if use_storey else "count_p_adj"
            fout.write(f"bin_start\tbin_end\tcount_p_value\t{adj_colname}\n")
        else:
            fout.write("bin_start\tbin_end\tcount_p_value\n")
        for i in range(bin_count):
            if i < fine_bins:
                start = i * fine_width
                end = start + fine_width
            else:
                start = 0.1 + (i - fine_bins) * coarse_width
                end = start + coarse_width
            if use_fdr:
                fout.write(f"{start:.4f}\t{end:.4f}\t{hist_counts_p[i]}\t{hist_counts_q[i]}\n")
            else:
                fout.write(f"{start:.4f}\t{end:.4f}\t{hist_counts_p[i]}\n")

    with detailed_path.open("w") as fout:
        if use_fdr:
            fout.write(
                f"unitig\t{count1_name}\t{count2_name}\t"
                f"{prev1_name}\t{prev2_name}\tglobal_prev\t"
                f"delta_prev\todds_ratio\tbias\tp_value\t{adj_label}\n"
            )
        else:
            fout.write(
                f"unitig\t{count1_name}\t{count2_name}\t"
                f"{prev1_name}\t{prev2_name}\tglobal_prev\t"
                f"delta_prev\todds_ratio\tbias\tp_value\n"
            )
        for row in stats_rows_sorted:
            def fmt(x):
                if isinstance(x, float) and math.isnan(x):
                    return "NA"
                if isinstance(x, float):
                    return f"{x:.6g}"
                return str(x)
            bias_label = row["bias"]
            if bias_label == "trait1":
                bias_label = trait1
            elif bias_label == "trait2":
                bias_label = trait2
            if use_fdr:
                fout.write(
                    f"{row['unitig']}\t"
                    f"{row['count_trait1']}\t{row['count_trait2']}\t"
                    f"{fmt(row['prev_trait1'])}\t{fmt(row['prev_trait2'])}\t"
                    f"{fmt(row['global_prev'])}\t{fmt(row['delta_prev'])}\t"
                    f"{fmt(row['odds_ratio'])}\t{bias_label}\t"
                    f"{fmt(row['p_value'])}\t{fmt(row.get(adj_label, math.nan))}\n"
                )
            else:
                fout.write(
                    f"{row['unitig']}\t"
                    f"{row['count_trait1']}\t{row['count_trait2']}\t"
                    f"{fmt(row['prev_trait1'])}\t{fmt(row['prev_trait2'])}\t"
                    f"{fmt(row['global_prev'])}\t{fmt(row['delta_prev'])}\t"
                    f"{fmt(row['odds_ratio'])}\t{bias_label}\t"
                    f"{fmt(row['p_value'])}\n"
                )

    logger.log(f"Summary and histograms written to {summary_path.name}")
    logger.log(f"Detailed stats written to {detailed_path.name}")
    logger.log("PASS V Fisher finished")
    logger.close()


if __name__ == "__main__":
    main()
