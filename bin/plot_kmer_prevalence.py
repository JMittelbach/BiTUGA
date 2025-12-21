#!/usr/bin/env python3

import argparse
import os
import sys
from pathlib import Path
import re
import csv
import math

try:
    import matplotlib.ticker as mtick
    from matplotlib.patches import Rectangle
except ImportError:
    mtick = None
    Rectangle = None

from plot_colors import get_color_map
from plot_utils import read_k

PALETTE = ["#c0392b", "#2471a3", "#8e44ad", "#16a085", "#d35400", "#2c3e50"]

_SUPERSCRIPTS = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")

def format_count(n: int) -> str:
    """Format counts als kompakte Zeichenkette mit ×10ⁿ (wie im Funnel-Plot-Stil)."""
    if n == 0:
        return "0"
    exp = int(math.floor(math.log10(abs(n))))
    if exp < 4:
        return str(n)
    mantissa = n / (10 ** exp)
    exp_str = str(exp).translate(_SUPERSCRIPTS)
    return f"{mantissa:.2g}×10{exp_str}"


def load_matplotlib(config_dir: Path):
    """Initialisiert Matplotlib im Funnel-Stil."""
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
        print("[WARN] matplotlib not installed; skipping k-mer prevalence plots.", file=sys.stderr)
        return None


def format_superscript_exponent(axis, font_size: int):
    """Konvertiert 1e+07 in ×10⁷ im Offset-Label einer Achse."""
    txt = axis.get_offset_text().get_text()
    if txt and ("e" in txt or "E" in txt):
        m = re.search(r"[eE]([+-]?\d+)", txt)
        if m:
            exp = m.group(1)
            axis.get_offset_text().set_text(f"×10{exp.translate(_SUPERSCRIPTS)}")
            axis.get_offset_text().set_fontsize(font_size)


def apply_sci_formatter(axis):
    """ScalarFormatter + MathText, damit das Offset-Label als ×10^n gesetzt werden kann."""
    fmt = mtick.ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((-4, 4))
    axis.set_major_formatter(fmt)


def parse_group_prevalence(path: Path):
    """
    Parsed die group_prevalence_histogram.txt und liefert zwei Dicts:
    - prevalence_data: {section: [(prevalence, n_kmers), ...]}
    - summaries: {key: value} aus den summary-Blöcken
    """
    prevalence_data = {}
    summaries = {}
    current = None
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("# "):
                if line.startswith("# "):
                    current = line[2:].strip()
                continue
            if current is None:
                continue
            parts = line.split("\t")
            if len(parts) == 2 and parts[0] and parts[1]:
                try:
                    x = int(parts[0])
                    y = int(parts[1])
                    prevalence_data.setdefault(current, []).append((x, y))
                    continue
                except ValueError:
                    pass
            if len(parts) == 2:
                key, val = parts
                try:
                    summaries[key] = int(val)
                except ValueError:
                    summaries[key] = val
    return prevalence_data, summaries


def write_hist_tsv(out_path: Path, trait_sections, prevalence_data: dict):
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["trait", "prevalence", "n_kmers"])
        for trait in trait_sections:
            for p, n in sorted(prevalence_data.get(trait, [])):
                writer.writerow([trait, p, n])


def write_diff_tsv(out_path: Path, diff_rows):
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["prevalence_diff", "n_kmers"])
        for p, n in diff_rows:
            writer.writerow([p, n])


def main():
    ap = argparse.ArgumentParser(description="Plot k-mer prevalence histograms/CDFs per trait.")
    ap.add_argument("--run-root", default=".", help="Run directory (default: .)")
    ap.add_argument(
        "--stats",
        help="Path to group_prevalence_histogram.txt "
        "(default: <run-root>/results/statistics/group_prevalence_histogram.txt)",
    )
    ap.add_argument(
        "--out-prefix",
        help="Output prefix (default: <run-root>/results/statistics/plots/kmer_stats/kmer_group_prevalence)",
    )
    ap.add_argument(
        "--palette",
        choices=["default", "mint"],
        default="default",
        help="Color palette for traits (default funnel colors; mint = soft blue/green tones).",
    )
    ap.add_argument(
        "--linear-y",
        action="store_true",
        help="Use linear y-scale for the histogram (default: log10(count)).",
    )
    args = ap.parse_args()

    root = Path(args.run_root).resolve()
    stats_path = (
        Path(args.stats)
        if args.stats
        else root / "results/statistics/group_prevalence_histogram.txt"
    )
    if not stats_path.is_file():
        sys.exit(f"[ERROR] Stats file not found: {stats_path}")

    out_prefix = (
        Path(args.out_prefix)
        if args.out_prefix
        else root / "results/statistics/plots/kmer_stats/kmer_group_prevalence"
    )
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    plt = load_matplotlib(out_prefix.parent / "mplconfig")
    if plt is None or mtick is None or Rectangle is None:
        sys.exit(0)
    prevalence_data, summaries = parse_group_prevalence(stats_path)

    diff_key = None
    for key in prevalence_data:
        if key.startswith("diff_abs("):
            diff_key = key
            break
    trait_sections = [k for k in prevalence_data.keys() if k != diff_key]
    if not trait_sections:
        sys.exit(f"[ERROR] No trait sections found in {stats_path}")
    trait_sections = sorted(trait_sections)
    has_diff = diff_key is not None

    font_size = 15
    legend_font_size = 14
    colors = get_color_map(trait_sections, palette=args.palette)
    k_val = read_k(root)
    k_label = f"{k_val}-mer" if k_val else r"$k$-mer"

    max_prev = 0
    fig_hist_width = 7.2
    width = 0.38
    all_counts = []
    zero_counts_present = False
    for i, trait in enumerate(trait_sections):
        rows = sorted(prevalence_data.get(trait, []))
        if not rows:
            continue
        prev, counts = zip(*rows)
        all_counts.extend(counts)
        if any(c == 0 for c in counts):
            zero_counts_present = True
        max_prev = max(max_prev, max(prev))
    fig_hist_width = 7.2
    if max_prev >= 15:
        width = 0.32
    if max_prev >= 30:
        width = 0.24

    fig_hist, ax_hist = plt.subplots(1, 1, figsize=(fig_hist_width, 4.6))
    for i, trait in enumerate(trait_sections):
        rows = sorted(prevalence_data.get(trait, []))
        if not rows:
            continue
        prev, counts = zip(*rows)
        offsets = [p + (i - 0.5) * width for p in prev]
        ax_hist.bar(
            offsets,
            counts,
            width=width,
            color=colors[trait],
            alpha=0.8,
            label=trait,
            edgecolor="none",
        )

    ax_hist.set_xlabel("Prevalence per trait group (samples)", fontsize=font_size)
    ax_hist.set_ylabel(fr"Unique {k_label}s", fontsize=font_size)
    if any(c > 0 for c in all_counts):
        ax_hist.set_yscale("log")
    all_ticks = list(range(0, max_prev + 1))
    ax_hist.set_xticks(all_ticks)
    if max_prev <= 10:
        ax_hist.set_xticklabels(all_ticks)
    elif max_prev <= 24:
        step = 2
        ax_hist.set_xticklabels([t if t % step == 0 else "" for t in all_ticks])
    else:
        step = 5
        ax_hist.set_xticklabels([t if t % step == 0 else "" for t in all_ticks])
    ax_hist.tick_params(axis="both", labelsize=font_size)
    legend_handles = [
        Rectangle((0, 0), 0.4, 0.4, color=colors[t], label=f"{t}-biased", alpha=0.8)
        for t in trait_sections
    ]
    ax_hist.legend(
        handles=legend_handles,
        frameon=False,
        fontsize=legend_font_size,
        handlelength=0.8,
        handleheight=0.8,
        loc="upper right",
    )
    run_label = root.name.replace("_", " ")
    ax_hist.set_title(
        f"Trait-biased unique $k$-mers across samples\n({run_label})",
        fontsize=font_size + 1,
        pad=18,
    )
    fig_hist.tight_layout()

    if has_diff:
        diff_rows = sorted(prevalence_data.get(diff_key, []))
        trait_a, trait_b = trait_sections[0], trait_sections[1] if len(trait_sections) > 1 else ("trait1", "trait2")
        m = re.search(r"diff_abs\(\|\s*([^-\|]+?)\s*-\s*([^-\|]+?)\s*\|\)", diff_key or "")
        if m:
            trait_a, trait_b = m.group(1).strip(), m.group(2).strip()
        fig_diff_h, ax_diff_h = plt.subplots(1, 1, figsize=(7.2, 4.0))
        if diff_rows:
            prev, counts = zip(*diff_rows)
            widths = [1.0] * len(prev)
            bars = ax_diff_h.bar(
                prev,
                counts,
                width=widths,
                color="#888888",
                edgecolor="none",
                alpha=0.9,
                align="center",
            )
            min_p, max_p = min(prev), max(prev)
            ticks = list(range(min_p, max_p + 1))
            ax_diff_h.set_xticks(ticks)
            if max_p <= 10:
                ax_diff_h.set_xticklabels(ticks)
            elif max_p <= 24:
                step = 2
                ax_diff_h.set_xticklabels([t if t % step == 0 else "" for t in ticks])
            else:
                step = 5
                ax_diff_h.set_xticklabels([t if t % step == 0 else "" for t in ticks])
            ax_diff_h.set_xlim(min_p - 1, max_p + 1)
            if any(c > 0 for c in counts):
                ax_diff_h.set_yscale("log")
        ax_diff_h.set_xlabel(
            f"Prevalence difference |{trait_a} - {trait_b}| (samples)", fontsize=font_size
        )
        ax_diff_h.set_ylabel(fr"Unique {k_label}s", fontsize=font_size)
        ax_diff_h.tick_params(axis="both", labelsize=font_size)
        ax_diff_h.set_title(
            f"Difference in trait-biased unique $k$-mers between groups\n({run_label})",
            fontsize=font_size + 1,
            pad=18,
        )
        fig_diff_h.tight_layout()

    hist_prefix = out_prefix.with_name(out_prefix.name + "_hist")
    fig_hist.savefig(hist_prefix.with_suffix(".pdf"), dpi=600)
    fig_hist.savefig(hist_prefix.with_suffix(".png"), dpi=600)
    if has_diff:
        diff_hist_prefix = out_prefix.with_name(out_prefix.name + "_diff_hist")
        fig_diff_h.savefig(diff_hist_prefix.with_suffix(".pdf"), dpi=600)
        fig_diff_h.savefig(diff_hist_prefix.with_suffix(".png"), dpi=600)
        plt.close(fig_diff_h)
    plt.close(fig_hist)

    hist_tsv = out_prefix.with_name(out_prefix.name + "_hist.tsv")
    write_hist_tsv(hist_tsv, trait_sections, prevalence_data)
    if has_diff:
        diff_hist_prefix = out_prefix.with_name(out_prefix.name + "_diff_hist")
        diff_tsv = diff_hist_prefix.with_suffix(".tsv")
        write_diff_tsv(diff_tsv, diff_rows)
    else:
        pass


if __name__ == "__main__":
    main()
