#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path
import math
import re
from plot_colors import get_color_map
from plot_utils import read_k

def log_warn(msg: str):
    print(f"[WARN] {msg}", file=sys.stderr)

def load_matplotlib(config_dir: Path):
    cfg_env = os.environ.get("MPLCONFIGDIR")
    cfg = Path(cfg_env) if cfg_env else Path(config_dir)
    try:
        cfg.mkdir(parents=True, exist_ok=True)
    except Exception:
        log_warn(f"Cannot create MPLCONFIGDIR at {cfg}; skipping unitig length plots.")
        return None
    os.environ["MPLCONFIGDIR"] = str(cfg)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
        from matplotlib.lines import Line2D
        matplotlib.rcParams["font.family"] = "sans-serif"
        matplotlib.rcParams["font.sans-serif"] = ["DejaVu Sans"]
        matplotlib.rcParams["pdf.fonttype"] = 42
        matplotlib.rcParams["ps.fonttype"] = 42
        return plt, Rectangle, Line2D
    except Exception as exc:
        log_warn(f"matplotlib not usable ({exc}); skipping unitig length plots.")
        return None

def parse_histogram(path: Path, header_marker: str):
    """Parse a length histogram block beginning with header_marker."""
    if not path.is_file():
        return []
    rows = []
    found = False
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                if found:
                    break
                continue
            if line.startswith("#") and header_marker in line:
                found = True
                continue
            if found:
                parts = line.split("\t")
                if len(parts) != 2:
                    continue
                try:
                    L = int(parts[0])
                    c = int(parts[1])
                except ValueError:
                    continue
                rows.append((L, c))
    return rows

def write_tsv(rows, out_path: Path, label: str):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        fh.write(f"{label}_len\tcount\n")
        for L, c in rows:
            fh.write(f"{L}\t{c}\n")

def lengths_from_fasta(fasta_path: Path):
    """
    Return (length histogram rows, id->length) for a FASTA/FASTA.GZ.
    Tries .fa/.fasta/.fa.gz/.fasta.gz fallbacks if the given path is missing.
    """
    candidates = [fasta_path]
    if fasta_path.suffix == ".fa":
        candidates.append(fasta_path.with_suffix(".fasta"))
        candidates.append(fasta_path.with_suffix(".fa.gz"))
        candidates.append(fasta_path.with_suffix(".fasta.gz"))
    elif fasta_path.suffix == ".fasta":
        candidates.append(fasta_path.with_suffix(".fa"))
        candidates.append(fasta_path.with_suffix(".fasta.gz"))
        candidates.append(fasta_path.with_suffix(".fa.gz"))
    else:
        candidates.append(fasta_path.with_suffix(fasta_path.suffix + ".gz"))

    existing = next((p for p in candidates if p.is_file()), None)
    if existing is None:
        log_warn(f"significant unitigs FASTA not found (tried {', '.join(str(p) for p in candidates)}); skipping significant plots.")
        return [], {}

    counts = {}
    id_to_len = {}
    current_id = None
    seq_len = 0

    import gzip
    opener = gzip.open if existing.suffix.endswith(".gz") else open
    with opener(existing, "rt", encoding="utf-8", errors="ignore") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_len > 0 and current_id is not None:
                    counts[seq_len] = counts.get(seq_len, 0) + 1
                    id_to_len[current_id] = seq_len
                current_id = line[1:].split()[0]
                seq_len = 0
            else:
                seq_len += len(line)
        if seq_len > 0 and current_id is not None:
            counts[seq_len] = counts.get(seq_len, 0) + 1
            id_to_len[current_id] = seq_len
    return sorted(counts.items()), id_to_len

def count_traits(metadata_path: Path):
    """Return (trait_counts dict, trait_column_name) or (None, None) on failure."""
    if not metadata_path.is_file():
        log_warn(f"metadata not found at {metadata_path}; skipping prevalence scatter")
        return None, None
    try:
        with metadata_path.open() as fh:
            header = fh.readline().strip().split("\t")
            if not header or len(header) < 2:
                return None, None
            trait_col_idx = None
            trait_keys = ["trait", "phenotype", "sex", "group"]
            header_lower = [h.lower() for h in header]
            for key in trait_keys:
                if key in header_lower:
                    trait_col_idx = header_lower.index(key)
                    break
            if trait_col_idx is None and len(header) >= 2:
                trait_col_idx = 1
            counts = {}
            for line in fh:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= trait_col_idx:
                    continue
                trait_val = parts[trait_col_idx]
                counts[trait_val] = counts.get(trait_val, 0) + 1
            return counts, header[trait_col_idx]
    except Exception as exc:
        log_warn(f"Failed to read metadata for trait counts ({exc}); skipping prevalence scatter")
        return None, None
    return None, None

def plot_hist(rows, title, ylabel, out_prefix: Path, plt, color="#555555", font_size=15):
    lengths, counts = zip(*sorted(rows))
    heights = [math.log10(c) if c > 0 else math.nan for c in counts]
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.bar(lengths, heights, width=0.9, color=color, edgecolor="none", alpha=0.9)
    ax.set_xlabel(title, fontsize=font_size)
    ax.set_ylabel(ylabel, fontsize=font_size)
    ax.tick_params(axis="both", labelsize=font_size)
    fig.tight_layout()
    fig.savefig(out_prefix.with_suffix(".pdf"), dpi=600)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=600)
    plt.close(fig)

def parse_length_classes(log_path: Path):
    """Return (short_max, medium_start, medium_end) if found in the unitig_matching log."""
    if not log_path.is_file():
        return (None, None, None)
    short_max = None
    medium_start = None
    medium_end = None
    pattern_short = re.compile(r"short\s+len\s*<=\s*(\d+)")
    pattern_medium = re.compile(r"medium\s+(\d+)\s*<=\s*len\s*<=\s*(\d+)")
    pattern_long = re.compile(r"long\s+len\s*>\s*(\d+)")
    with log_path.open() as fh:
        for line in fh:
            if "short" in line and "len" in line:
                m = pattern_short.search(line)
                if m:
                    short_max = int(m.group(1))
            if "medium" in line and "len" in line:
                m = pattern_medium.search(line)
                if m:
                    medium_start = int(m.group(1))
                    medium_end = int(m.group(2))
            if "long" in line and "len" in line and medium_end is None:
                m = pattern_long.search(line)
                if m and medium_start is not None:
                    # If only long threshold is given, infer medium_end from it
                    medium_end = int(m.group(1))
    return (short_max, medium_start, medium_end)

def plot_hist_with_classes(
    rows,
    xlabel,
    ylabel,
    out_prefix,
    plt,
    Rectangle,
    Line2D,
    short_max,
    medium_start,
    medium_end,
    font_size=15,
    pad=None,
    x_min=None,
    x_max=None,
    use_classes=True,
    title=None,
):
    lengths, counts = zip(*sorted(rows))
    if not any(c > 0 for c in counts):
        log_warn(f"No positive counts for {out_prefix.name}; skipping plot.")
        return
    colors = []
    col_small = "#c8cbcf"
    col_medium = "#9aa0a6"
    col_long = "#5b626a"
    if use_classes:
        for L in lengths:
            if short_max is not None and L <= short_max:
                colors.append(col_small)
            elif medium_start is not None and medium_end is not None:
                if L < medium_start:
                    colors.append(col_small)
                elif L <= medium_end:
                    colors.append(col_medium)
                else:
                    colors.append(col_long)
            else:
                colors.append(col_long)
    else:
        colors = ["#888888"] * len(lengths)
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.bar(lengths, counts, width=1.0, color=colors, edgecolor="none", alpha=0.95, align="center")
    ax.set_yscale("log")
    xmin_local = min(lengths)
    xmax_local = max(lengths)
    xmin = x_min if x_min is not None else xmin_local
    xmax = x_max if x_max is not None else xmax_local
    pad_eff = pad
    if pad_eff is None:
        pad_eff = max(1.0, 0.02 * (xmax_local - xmin_local))
    ax.set_xlim(xmin - pad_eff, xmax + pad_eff)
    ax.set_xlabel(xlabel, fontsize=font_size)
    ax.set_ylabel(ylabel, fontsize=font_size)
    if title:
        ax.set_title(title, fontsize=font_size + 1, pad=18)
    ax.tick_params(axis="both", labelsize=font_size)

    total_counts = sum(counts)
    mean_len = None
    if total_counts > 0:
        mean_len = sum(l * c for l, c in zip(lengths, counts)) / total_counts
        ax.axvline(mean_len, color="#333333", linestyle="--", linewidth=1.0)

    legend_handles = []
    labels = []
    if use_classes:
        if short_max is not None:
            legend_handles.append(Rectangle((0, 0), 0.4, 0.4, color=col_small, ec="none"))
            labels.append(f"short (≤{short_max} bp)")
        if medium_start is not None and medium_end is not None:
            legend_handles.append(Rectangle((0, 0), 0.4, 0.4, color=col_medium, ec="none"))
            labels.append(f"medium ({medium_start}–{medium_end} bp)")
        legend_handles.append(Rectangle((0, 0), 0.4, 0.4, color=col_long, ec="none"))
        long_label = "long"
        if medium_end is not None:
            long_label += f" (> {medium_end} bp)"
        elif short_max is not None:
            long_label += f" (> {short_max} bp)"
        labels.append(long_label)
    if mean_len is not None:
        legend_handles.append(Line2D([0], [0], color="#333333", linestyle="--", linewidth=1.0))
        labels.append("mean")
    ax.legend(legend_handles, labels, loc="upper right", frameon=False, fontsize=font_size - 1, handlelength=0.8, handleheight=0.8)

    fig.tight_layout()
    fig.savefig(out_prefix.with_suffix(".pdf"), dpi=600)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=600)
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser(description="Plot unitig length histograms (pre- and post-matching).")
    ap.add_argument("--run-root", default=".", help="Run directory (default: .)")
    ap.add_argument("--stats", help="Path to unitig_statistics.txt (default: <run-root>/results/statistics/unitig_statistics.txt)")
    ap.add_argument("--matches", help="Path to unitig_matches.txt (default: <run-root>/results/statistics/unitig_matches.txt)")
    ap.add_argument("--significant", help="Path to significant_unitigs.fa (default: <run-root>/results/significant_unitigs.fa)")
    ap.add_argument("--out-prefix", help="Output prefix (default: <run-root>/results/statistics/plots/unitigs/unitig_lengths)")
    ap.add_argument("--matching-log", help="Path to unitig_matching.log (default: <run-root>/logging/unitig_matching.log)")
    ap.add_argument("--metadata", help="Path to metadata.tsv (default: <run-root>/metadata.tsv)")
    ap.add_argument("--presence", help="Path to unitig_presence.tsv (default: <run-root>/results/unitig_presence.tsv)")
    args = ap.parse_args()

    root = Path(args.run_root).resolve()
    stats_path = Path(args.stats) if args.stats else root / "results/statistics/unitig_statistics.txt"
    matches_path = Path(args.matches) if args.matches else root / "results/statistics/unitig_matches.txt"
    signif_path = Path(args.significant) if args.significant else root / "results/significant_unitigs.fa"
    out_prefix = Path(args.out_prefix) if args.out_prefix else root / "results/statistics/plots/unitigs/unitig_lengths"
    cfg_dir = root / "results/statistics/plots/mplconfig"
    matching_log = Path(args.matching_log) if args.matching_log else root / "logging/unitig_matching.log"
    metadata_path = Path(args.metadata) if args.metadata else root / "metadata.tsv"
    presence_path = Path(args.presence) if args.presence else root / "results/unitig_presence.tsv"
    k_val = read_k(root)
    k_label = f"{k_val}-mer" if k_val else r"$k$-mer"

    mpl = load_matplotlib(cfg_dir)
    if mpl is None:
        return
    plt, Rectangle, Line2D = mpl

    all_rows = parse_histogram(stats_path, "length_histogram")
    present_rows = parse_histogram(matches_path, "present unitigs length histogram")
    signif_rows, signif_len_map = lengths_from_fasta(signif_path)

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(all_rows, out_prefix.with_name(out_prefix.name + "_candidates.tsv"), "length")
    write_tsv(present_rows, out_prefix.with_name(out_prefix.name + "_present.tsv"), "length")
    write_tsv(signif_rows, out_prefix.with_name(out_prefix.name + "_significant.tsv"), "length")

    short_max, medium_start, medium_end = parse_length_classes(matching_log)
    shared_pad = None
    combined_lengths = []
    if all_rows:
        combined_lengths.extend([L for L, _ in all_rows])
    if present_rows:
        combined_lengths.extend([L for L, _ in present_rows])
    if signif_rows:
        combined_lengths.extend([L for L, _ in signif_rows])
    if all_rows:
            plot_hist_with_classes(
                all_rows,
                xlabel="Unitig length (bp)",
                ylabel="Unitigs (candidates)",
                out_prefix=out_prefix.with_name(out_prefix.name + "_candidates"),
                plt=plt,
                Rectangle=Rectangle,
                Line2D=Line2D,
                short_max=short_max,
                medium_start=medium_start,
                medium_end=medium_end,
                title=f"Candidate unitig length distribution\n({root.name.replace('_', ' ')})",
            )
    else:
        log_warn("No histogram data found in unitig_statistics.txt; skipping all-unitigs plot.")

    if present_rows:
            plot_hist_with_classes(
                present_rows,
                xlabel="Unitig length (bp)",
                ylabel="Unitigs (matching)",
                out_prefix=out_prefix.with_name(out_prefix.name + "_present"),
                plt=plt,
                Rectangle=Rectangle,
                Line2D=Line2D,
                short_max=short_max,
                medium_start=medium_start,
                medium_end=medium_end,
                use_classes=False,
                title=f"Matching unitig length distribution\n({root.name.replace('_', ' ')})",
            )
    else:
        log_warn("No histogram data found in unitig_matches.txt; skipping present-unitigs plot.")

    if signif_rows:
            plot_hist_with_classes(
                signif_rows,
                xlabel="Unitig length (bp)",
                ylabel="Unitigs (significant)",
                out_prefix=out_prefix.with_name(out_prefix.name + "_significant"),
                plt=plt,
                Rectangle=Rectangle,
                Line2D=Line2D,
                short_max=short_max,
                medium_start=medium_start,
                medium_end=medium_end,
                use_classes=False,
                title=f"Significant unitig length distribution\n({root.name.replace('_', ' ')})",
            )
    else:
        log_warn("No significant_unitigs.fa found or empty; skipping significant-unitigs plot.")

    # silence noisy path logging

if __name__ == "__main__":
    main()
