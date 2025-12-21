#!/usr/bin/env python3

import argparse
import os
import sys
from pathlib import Path
import math
import re

try:
    import matplotlib.ticker as mtick
    from matplotlib.patches import Rectangle
except ImportError:
    mtick = None
    Rectangle = None

from plot_colors import get_color_map
from plot_utils import read_k

_SUPERSCRIPTS = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")

def format_superscript_exponent(axis, font_size):
    """
    Formatiert den Achsen-Offset-Text (z. B. 1e+07) in den Stil ×10ⁿ
    mit Superscripts, um den Stil des Funnel Plots zu matchen.
    """
    txt = axis.get_offset_text().get_text()
    if txt and ("e" in txt or "E" in txt):
        try:
            match = re.search(r"[Ee]([+-]?\d+)", txt)
            if match:
                exp = int(match.group(1))
                
                exp_str = str(exp).translate(_SUPERSCRIPTS)
                axis.get_offset_text().set_text(f"×10{exp_str}")
                axis.get_offset_text().set_fontsize(font_size) 
            else:
                 axis.get_offset_text().set_fontsize(font_size) 
        except Exception:
            pass 

def load_matplotlib(config_dir: Path):
    """Initialisiert die matplotlib-Konfiguration, um den Stil von funnel_plot.py zu matchen."""
    cfg = Path(os.environ.get("MPLCONFIGDIR", config_dir))
    cfg.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(cfg)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        matplotlib.rcParams["font.family"] = "sans-serif"
        matplotlib.rcParams["font.sans-serif"] = ["DejaVu Sans"]
        matplotlib.rcParams["pdf.fonttype"] = 42
        matplotlib.rcParams["ps.fonttype"] = 42
        return plt
    except ImportError:
        print("[WARN] matplotlib not installed; skipping k-mer stats plots.", file=sys.stderr)
        return None


def main():
    ap = argparse.ArgumentParser(description="Plot k-mer statistics per sample.")
    ap.add_argument("--run-root", default=".", help="Run directory (default: .)")
    ap.add_argument(
        "--stats",
        help="Path to kmer_stats_per_sample.tsv "
        "(default: <run-root>/results/statistics/kmer_stats_per_sample.tsv)",
    )
    ap.add_argument(
        "--out-prefix",
        help="Output prefix (default: <run-root>/results/statistics/plots/kmer_stats/kmer_stats)",
    )
    ap.add_argument(
        "--palette",
        choices=["default", "mint"],
        default="default",
        help="Color palette for traits (default funnel colors; mint = soft blue/green tones).",
    )
    args = ap.parse_args()

    root = Path(args.run_root).resolve()
    stats_path = Path(args.stats) if args.stats else root / "results/statistics/kmer_stats_per_sample.tsv"
    if not stats_path.is_file():
        sys.exit(f"[ERROR] Stats file not found: {stats_path}")

    out_prefix = (
        Path(args.out_prefix)
        if args.out_prefix
        else root / "results/statistics/plots/kmer_stats/kmer_stats"
    )
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    plt = load_matplotlib(out_prefix.parent / "mplconfig")
    if plt is None or mtick is None or Rectangle is None:
        sys.exit(0)

    import csv

    rows = []
    with stats_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    if not rows:
        sys.exit(f"[ERROR] No data in {stats_path}")

    required_cols = {"sample_id", "trait", "unique_counted_kmers", "total_kmers", "total_reads"}
    if not required_cols.issubset(rows[0].keys()):
        missing = required_cols - set(rows[0].keys())
        sys.exit(f"[ERROR] Missing columns in {stats_path}: {', '.join(sorted(missing))}")

    def to_int(val):
        try:
            return int(val)
        except Exception:
            try:
                return int(float(val))
            except Exception:
                return 0

    data = []
    for r in rows:
        data.append(
            {
                "sample_id": r.get("sample_id", ""),
                "trait": r.get("trait", ""),
                "unique_counted_kmers": to_int(r.get("unique_counted_kmers", 0)),
                "total_kmers": to_int(r.get("total_kmers", 0)),
                "total_reads": to_int(r.get("total_reads", 0)),
            }
        )

    trait_order = sorted({d["trait"] for d in data if d["trait"]})
    if not trait_order:
        sys.exit("[ERROR] No traits found in stats file.")

    font_size = 15
    legend_font_size = 14
    colors = get_color_map(trait_order, palette=args.palette)
    k_val = read_k(root)
    k_label = f"{k_val}-mer" if k_val else r"$k$-mer"

    trait_values = {t: [] for t in trait_order}
    for d in data:
        trait_values[d["trait"]].append(d["unique_counted_kmers"])

    def apply_sci_formatter(axis):
        """Nutze ScalarFormatter mit MathText, damit Offsets als ×10^n erscheinen."""
        fmt = mtick.ScalarFormatter(useMathText=True)
        fmt.set_powerlimits((-4, 4))
        axis.set_major_formatter(fmt)

    fig_box, ax_box = plt.subplots(1, 1, figsize=(7.2, 4.8))
    positions = [i * 0.7 for i in range(len(trait_order))]
    
    for pos, t in zip(positions, trait_order):
        vals_t = trait_values[t]
        if len(vals_t) >= 10:
            vparts = ax_box.violinplot(
                [vals_t],
                positions=[pos],
                showmeans=False,
                showmedians=False,
                showextrema=False,
                widths=0.6,
            )
            for body in vparts["bodies"]:
                body.set_facecolor(colors[t])
                body.set_edgecolor(colors[t])
                body.set_alpha(0.35)
    bplot = ax_box.boxplot(
        [trait_values[t] for t in trait_order],
        positions=positions,
        widths=0.4, 
        patch_artist=True,
        showfliers=False,
    )
    for patch, t in zip(bplot["boxes"], trait_order):
        patch.set_facecolor("none")
        patch.set_edgecolor("black")
        patch.set_linewidth(1.5)
    for element in ["whiskers", "caps", "medians"]:
        for line in bplot[element]:
            line.set_color("black")
            line.set_linewidth(1.2)
            
    jitter = 0.05
    for pos, t in zip(positions, trait_order):
        ys = trait_values[t]
        xs = [pos + (jitter * (2 * (i % 2) - 1)) * 0.25 for i in range(len(ys))]
        ax_box.scatter(xs, ys, s=30, alpha=0.7, color=colors[t], edgecolors="none")

    ax_box.set_xticks(positions)
    ax_box.set_xticklabels(trait_order, rotation=0, fontsize=font_size)
    ax_box.set_xlabel("Trait", fontsize=font_size)
    ax_box.set_ylabel(fr"Unique {k_label}s", fontsize=font_size)
    run_label = root.name.replace("_", " ")
    ax_box.set_title(
        f"Distribution of unique $k$-mers per trait\n({run_label})",
        fontsize=font_size + 1,
        pad=26,
    )
    ax_box.tick_params(axis="both", labelsize=font_size)
    apply_sci_formatter(ax_box.yaxis)
    format_superscript_exponent(ax_box.yaxis, font_size)

    fig_scat, ax_scat = plt.subplots(1, 1, figsize=(7.2, 4.8))
    for t in trait_order:
        xs = [d["total_reads"] for d in data if d["trait"] == t]
        ys = [d["unique_counted_kmers"] for d in data if d["trait"] == t]
        ax_scat.scatter(xs, ys, s=40, alpha=0.8, color=colors[t], edgecolors="none")
        
    ax_scat.set_xlabel("Total reads", fontsize=font_size)
    ax_scat.set_ylabel(fr"Unique {k_label}s", fontsize=font_size)
    ax_scat.tick_params(axis="both", labelsize=font_size)
    apply_sci_formatter(ax_scat.xaxis)
    apply_sci_formatter(ax_scat.yaxis)
    format_superscript_exponent(ax_scat.xaxis, font_size)
    format_superscript_exponent(ax_scat.yaxis, font_size)
    ax_scat.set_title(
        f"Unique $k$-mers vs total reads\n({run_label})",
        fontsize=font_size + 1,
        pad=18,
    )
    
    legend_handles = [
        Rectangle((0, 0), 0.4, 0.4, color=colors[t], label=t, alpha=0.8)
        for t in trait_order
    ]
    ax_scat.legend(
        handles=legend_handles, 
        frameon=False, 
        fontsize=legend_font_size, 
        handlelength=0.8, 
        handleheight=0.8,
        loc="upper left"
    )

    fig_box.tight_layout()
    fig_scat.tight_layout()
    box_prefix = out_prefix.with_name(out_prefix.name + "_box")
    scat_prefix = out_prefix.with_name(out_prefix.name + "_scatter")
    fig_box.savefig(box_prefix.with_suffix(".pdf"), dpi=600)
    fig_box.savefig(box_prefix.with_suffix(".png"), dpi=600)
    fig_scat.savefig(scat_prefix.with_suffix(".pdf"), dpi=600)
    fig_scat.savefig(scat_prefix.with_suffix(".png"), dpi=600)
    plt.close(fig_box)
    plt.close(fig_scat)

    out_tsv = out_prefix.with_suffix(".tsv")
    with out_tsv.open("w") as fh:
        fh.write("sample_id\ttrait\tunique_counted_kmers\ttotal_kmers\ttotal_reads\n")
        for d in data:
            fh.write(
                f"{d['sample_id']}\t{d['trait']}\t{d['unique_counted_kmers']}\t{d['total_kmers']}\t{d['total_reads']}\n"
            )

    def quantiles(vals):
        vals_sorted = sorted(vals)
        n = len(vals_sorted)
        if n == 0:
            return (math.nan, math.nan, math.nan)
        def percentile(p):
            k = (n - 1) * p
            f = math.floor(k)
            c = math.ceil(k)
            if f == c:
                return vals_sorted[int(k)]
            d0 = vals_sorted[int(f)] * (c - k)
            d1 = vals_sorted[int(c)] * (k - f)
            return d0 + d1
        q1 = percentile(0.25)
        q2 = percentile(0.50)
        q3 = percentile(0.75)
        return q1, q2, q3

    box_points_tsv = box_prefix.with_name(box_prefix.name + "_points.tsv")
    box_summary_tsv = box_prefix.with_name(box_prefix.name + "_summary.tsv")

    with box_points_tsv.open("w") as fh_pts, box_summary_tsv.open("w") as fh_sum:
        fh_pts.write("sample_id\ttrait\tunique_counted_kmers\tis_outlier\n")
        fh_sum.write("trait\tn\tmean\tmedian\tq1\tq3\tiqr\twhisker_low\twhisker_high\tn_outliers\n")
        for t in trait_order:
            vals = trait_values[t]
            n = len(vals)
            if n == 0:
                continue
            q1, med, q3 = quantiles(vals)
            iqr = q3 - q1
            whisker_low = min(v for v in vals if v >= q1 - 1.5 * iqr) if n > 0 else math.nan
            whisker_high = max(v for v in vals if v <= q3 + 1.5 * iqr) if n > 0 else math.nan
            outliers = [v for v in vals if v < whisker_low or v > whisker_high]
            mean_val = sum(vals) / n
            for d in [d for d in data if d["trait"] == t]:
                v = d["unique_counted_kmers"]
                is_out = (v < whisker_low) or (v > whisker_high)
                fh_pts.write(f"{d['sample_id']}\t{t}\t{v}\t{str(is_out).lower()}\n")
            fh_sum.write(
                f"{t}\t{n}\t{mean_val:.6g}\t{med:.6g}\t{q1:.6g}\t{q3:.6g}\t{iqr:.6g}\t{whisker_low:.6g}\t{whisker_high:.6g}\t{len(outliers)}\n"
            )



if __name__ == "__main__":
    main()
