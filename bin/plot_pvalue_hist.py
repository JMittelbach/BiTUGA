#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path
import math
import csv

from plot_utils import read_k

def log_warn(msg: str):
    print(f"[WARN] {msg}", file=sys.stderr)


def load_matplotlib(cfg_dir: Path):
    cfg = Path(os.environ.get("MPLCONFIGDIR", cfg_dir))
    try:
        cfg.mkdir(parents=True, exist_ok=True)
    except Exception:
        log_warn(f"Cannot create MPLCONFIGDIR at {cfg}; skipping p-value plots.")
        return None
    os.environ["MPLCONFIGDIR"] = str(cfg)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        matplotlib.rcParams["font.family"] = "sans-serif"
        matplotlib.rcParams["font.sans-serif"] = ["DejaVu Sans"]
        matplotlib.rcParams["pdf.fonttype"] = 42
        matplotlib.rcParams["ps.fonttype"] = 42
        return plt, Line2D
    except Exception as exc:
        log_warn(f"matplotlib not usable ({exc}); skipping p-value plots.")
        return None


def read_values(details_path: Path):
    p_vals = []
    adj_vals = []
    adj_label = None
    if not details_path.is_file():
        log_warn(f"Details file not found: {details_path}")
        return p_vals, adj_vals, adj_label
    with details_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            return p_vals, adj_vals, adj_label
        fields_lower = {f.lower(): f for f in reader.fieldnames}
        p_key = fields_lower.get("p_value") or fields_lower.get("p")
        if "q_value" in fields_lower:
            adj_key = fields_lower["q_value"]
            adj_label = "q_value"
        elif "p_adj" in fields_lower:
            adj_key = fields_lower["p_adj"]
            adj_label = "p_adj"
        else:
            adj_key = None
        if p_key is None:
            log_warn("No p_value column found in details file.")
            return p_vals, adj_vals, adj_label
        for row in reader:
            try:
                p = float(row[p_key])
                if p > 0:
                    p_vals.append(p)
            except Exception:
                pass
            if adj_key:
                try:
                    q = float(row[adj_key])
                    if q > 0:
                        adj_vals.append(q)
                except Exception:
                    pass
    return p_vals, adj_vals, adj_label


def make_hist(values, bins):
    if not values:
        return []
    bins_sorted = [i / bins for i in range(bins + 1)]
    counts = [0] * bins
    for v in values:
        if v < 0:
            continue
        if v > 1:
            v = 1.0
        idx = min(bins - 1, int(v * bins))
        counts[idx] += 1
    hist = []
    for i in range(bins):
        hist.append((bins_sorted[i], bins_sorted[i + 1], counts[i]))
    return hist


def plot_hist(
    hist,
    out_prefix: Path,
    label: str,
    plt,
    Line2D,
    font_size=14,
    threshold=None,
    threshold_label=None,
    title=None,
):
    if not hist:
        log_warn(f"No data for {label} histogram; skipping.")
        return
    fig, ax = plt.subplots(figsize=(7.0, 4.0))

    starts = [x[0] for x in hist]
    widths = [x[1] - x[0] for x in hist]
    counts = [x[2] for x in hist]
    ax.bar(starts, counts, width=widths, align="edge", color="#888888", edgecolor="none", alpha=0.9)

    ax.set_xlabel(label, fontsize=font_size)
    ax.set_ylabel("Tested unitigs", fontsize=font_size)
    ax.tick_params(axis="both", labelsize=font_size - 1)
    ax.set_xlim(0, 1)
    if title:
        ax.set_title(title, fontsize=font_size + 1, pad=16)
    if any(c > 0 for c in counts):
        ax.set_yscale("log")
    legend_handles = []
    legend_labels = []
    if threshold is not None and threshold > 0:
        ax.axvline(threshold, color="#333333", linestyle="--", linewidth=1.0)
        legend_handles.append(Line2D([0], [0], color="#333333", linestyle="--", linewidth=1.0))
        legend_labels.append(threshold_label or f"{threshold}")
    if legend_handles:
        ax.legend(
            handles=legend_handles,
            labels=legend_labels,
            frameon=False,
            fontsize=font_size - 1,
            handlelength=0.8,
            handleheight=0.8,
            loc="upper right",
        )
    fig.tight_layout()
    fig.savefig(out_prefix.with_suffix(".pdf"), dpi=600)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=600)
    plt.close(fig)


def compute_median_from_hist(hist):
    total = sum(c for _, _, c in hist)
    if total == 0:
        return 0.5
    cum = 0
    for start, end, count in hist:
        cum += count
        if cum >= total / 2:
            return (start + end) / 2
    return 0.5


def write_hist_tsv(hist, out_path: Path, label: str):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        fh.write("bin_start\tbin_end\tcount\n")
        for start, end, cnt in hist:
            fh.write(f"{start:.4f}\t{end:.4f}\t{cnt}\n")


def read_stats(stats_path: Path):
    """Return p_threshold, alpha if available."""
    if not stats_path.is_file():
        return None, None
    p_thr = None
    alpha = None
    with stats_path.open() as fh:
        for line in fh:
            if not line.strip() or "\t" not in line:
                continue
            key, val = line.rstrip("\n").split("\t", 1)
            if key == "p_threshold":
                try:
                    p_thr = float(val)
                except Exception:
                    pass
            if key == "alpha":
                try:
                    alpha = float(val)
                except Exception:
                    pass
    return p_thr, alpha


def main():
    ap = argparse.ArgumentParser(description="Plot p-value / adjusted p-value histograms.")
    ap.add_argument("--run-root", default=".", help="Run directory (default: .)")
    ap.add_argument("--details", help="unitig_association_details.tsv (default: <run-root>/results/statistics/unitig_association_details.tsv)")
    ap.add_argument("--out-prefix", help="Output prefix (default: <run-root>/results/statistics/plots/pvalues/pvalues)")
    ap.add_argument("--bins", type=int, default=50, help="Number of bins (default: 50)")
    ap.add_argument("--stats", help="unitig_association_stats.txt (default: <run-root>/results/statistics/unitig_association_stats.txt)")
    args = ap.parse_args()

    root = Path(args.run_root).resolve()
    run_label = root.name.replace("_", " ")
    details_path = Path(args.details) if args.details else root / "results/statistics/unitig_association_details.tsv"
    out_prefix = Path(args.out_prefix) if args.out_prefix else root / "results/statistics/plots/pvalues/pvalues"
    stats_path = Path(args.stats) if args.stats else root / "results/statistics/unitig_association_stats.txt"
    cfg_dir = root / "results/statistics/plots/mplconfig"
    k_val = read_k(root)
    k_label = f"{k_val}-mer" if k_val else r"$k$-mer"

    mpl = load_matplotlib(cfg_dir)
    if mpl is None:
        return
    plt, Line2D = mpl

    p_vals, adj_vals, adj_label = read_values(details_path)

    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    p_thr, alpha = read_stats(stats_path)

    max_bins = 200
    bins_p = min(max_bins, len(p_vals)) if len(p_vals) > 0 else min(max_bins, args.bins)
    hist_p = make_hist(p_vals, bins_p)
    write_hist_tsv(hist_p, out_prefix.with_name(out_prefix.name + "_p.tsv"), "p_value")
    plot_hist(
        hist_p,
        out_prefix.with_name(out_prefix.name + "_p"),
        "p-value",
        plt,
        Line2D,
        threshold=p_thr,
        threshold_label=f"p-threshold={p_thr}" if p_thr else None,
        title=f"Distribution of p-values for tested unitigs\n({run_label})",
    )

    if adj_vals:
        bins_adj = min(max_bins, len(adj_vals)) if len(adj_vals) > 0 else min(max_bins, args.bins)
        hist_adj = make_hist(adj_vals, bins_adj)
        write_hist_tsv(hist_adj, out_prefix.with_name(out_prefix.name + "_adj.tsv"), adj_label or "p_adj")
        if adj_label == "q_value":
            x_label = "q-value"
            title = f"Distribution of q-values for tested unitigs\n({run_label})"
        else:
            x_label = "adjusted p-value"
            title = f"Distribution of adjusted p-values for tested unitigs\n({run_label})"
        plot_hist(
            hist_adj,
            out_prefix.with_name(out_prefix.name + "_adj"),
            x_label,
            plt,
            Line2D,
            threshold=alpha,
            threshold_label=f"{adj_label or 'alpha'}={alpha}" if alpha else None,
            title=title,
        )


if __name__ == "__main__":
    main()
