#!/usr/bin/env python3

import argparse
import gzip
import math
import sys
import re
from pathlib import Path
import os
import traceback
from plot_colors import get_color_map
from plot_utils import read_k

def load_matplotlib(config_dir: Path):
    """Load matplotlib safely; return (plt, Rectangle) or (None, None) if missing."""
    global plt, Rectangle
    cfg = Path(os.environ.get("MPLCONFIGDIR", config_dir))
    cfg.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(cfg)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
        matplotlib.rcParams["font.family"] = "sans-serif"
        matplotlib.rcParams["font.sans-serif"] = ["DejaVu Sans"]
        matplotlib.rcParams["pdf.fonttype"] = 42
        matplotlib.rcParams["ps.fonttype"] = 42
        return plt, Rectangle
    except ImportError:
        print("[WARN] matplotlib not installed; skipping pipeline funnel plots.", file=sys.stderr)
        return None, None

def open_maybe_gz(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def count_fasta_records(path):
    if not path or not path.is_file():
        return 0
    with open_maybe_gz(path, "rt") as fh:
        return sum(1 for line in fh if line.startswith(">"))

def parse_kmer_union(stats_dir):
    def read_union(path):
        if not path.is_file():
            return None
        with path.open() as fh:
            for line in fh:
                if line.startswith("union_total"):
                    try:
                        return int(line.strip().split("\t")[1])
                    except Exception:
                        return None
        return None

    for cand in ["kmer_candidates.txt", "group_prevalence_histogram.txt", "kmer_totals.tsv"]:
        val = read_union(stats_dir / cand)
        if val is not None:
            return val
    return 0

def parse_unitig_stats(stats_dir):
    path = stats_dir / "unitig_statistics.txt"
    out = {}
    if not path.is_file():
        return out
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            if "\t" not in line:
                continue
            key, value = line.strip().split("\t", 1)
            try:
                out[key] = int(value)
            except ValueError:
                try:
                    out[key] = float(value)
                except ValueError:
                    out[key] = value
    return out

def count_presence(path):
    if not path.is_file():
        return 0
    with path.open() as fh:
        next(fh, None)
        return sum(1 for _ in fh)

def parse_fisher_summary(stats_dir):
    summary = stats_dir / "unitig_association_stats.txt"
    if not summary.is_file():
        return "unknown", None, {}, 0
    mode = "unknown"
    sig_total = None
    bias = {}
    bias_tie = 0
    with summary.open() as fh:
        for line in fh:
            if line.startswith("multiple_testing_correction"):
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    mode = parts[1].lower()
            elif line.startswith("significant_unitigs"):
                try:
                    sig_total = int(line.strip().split("\t")[1])
                except Exception:
                    sig_total = None
            elif line.startswith("significant_bias_"):
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    trait = parts[0].replace("significant_bias_", "")
                    try:
                        bias[trait] = int(parts[1])
                    except Exception:
                        bias[trait] = 0
            elif line.startswith("significant_bias_tie"):
                try:
                    bias_tie = int(line.strip().split("\t")[1])
                except Exception:
                    bias_tie = 0
    return mode, sig_total, bias, bias_tie

def parse_unitig_match_bias(stats_dir):
    path = stats_dir / "unitig_matches.txt"
    if not path.is_file():
        return {}, 0
    bias = {}
    bias_t = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith("TRAIT_BIAS"):
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    trait = parts[1]
                    vals = dict(p.split("=", 1) for p in parts[2:])
                    cnt = 0
                    try:
                        cnt = int(vals.get("specific_only", 0)) + int(vals.get("biased_gt_other", 0))
                    except Exception:
                        cnt = 0
                    bias[trait] = cnt
            elif line.startswith("abs_diff_hist"):
                for hist_line in fh:
                    hist_line = hist_line.strip()
                    if not hist_line:
                        continue
                    if "\t" not in hist_line:
                        break
                    try:
                        diff, n = hist_line.split("\t", 1)
                        if diff == "0":
                            bias_t = int(n)
                            break
                    except Exception:
                        continue
    return bias, bias_t

def parse_candidate_bias(stats_dir):
    def read_bias(path):
        bias_local = {}
        bias_t_local = 0
        if not path.is_file():
            return bias_local, bias_t_local
        with path.open() as fh:
            for line in fh:
                if line.startswith("post_bias_"):
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        trait = parts[0].replace("post_bias_", "")
                        try:
                            bias_local[trait] = int(parts[1])
                        except Exception:
                            bias_local[trait] = 0
                elif line.startswith("post_ties"):
                    try:
                        bias_t_local = int(line.strip().split("\t")[1])
                    except Exception:
                        bias_t_local = 0
        return bias_local, bias_t_local

    bias, bias_t = read_bias(stats_dir / "kmer_candidates.txt")
    if bias or bias_t:
        return bias, bias_t
    return read_bias(stats_dir / "group_prevalence_histogram.txt")

def parse_pre_bias(stats_dir):
    def read_pre(path):
        bias_local = {}
        bias_t_local = 0
        if not path.is_file():
            return bias_local, bias_t_local
        with path.open() as fh:
            for line in fh:
                if line.startswith("pre_bias_"):
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        trait = parts[0].replace("pre_bias_", "")
                        try:
                            bias_local[trait] = int(parts[1])
                        except Exception:
                            bias_local[trait] = 0
                elif line.startswith("pre_ties"):
                    try:
                        bias_t_local = int(line.strip().split("\t")[1])
                    except Exception:
                        bias_t_local = 0
        return bias_local, bias_t_local

    bias, bias_t = read_pre(stats_dir / "kmer_candidates.txt")
    if bias or bias_t:
        return bias, bias_t
    return read_pre(stats_dir / "group_prevalence_histogram.txt")

def parse_fisher(detailed_path, alpha, p_threshold, fdr_enabled):
    if not detailed_path.is_file():
        return 0, 0, 0, False
    with detailed_path.open() as fh:
        header = fh.readline().strip().split("\t")
        col = {name: idx for idx, name in enumerate(header)}
        idx_p = col.get("p") if col.get("p") is not None else col.get("p_value")
        idx_q = col.get("p_adj") if col.get("p_adj") is not None else col.get("q_value")
        tested = sig_p = sig_fdr = 0
        fdr_used = False
        for line in fh:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            tested += 1
            try:
                p = float(parts[idx_p]) if idx_p is not None else math.inf
            except Exception:
                p = math.inf
            try:
                q = float(parts[idx_q]) if idx_q is not None else math.inf
            except Exception:
                q = math.inf
            if p <= p_threshold:
                sig_p += 1
            if fdr_enabled and math.isfinite(q) and q <= alpha:
                sig_fdr += 1
                fdr_used = True
    return tested, sig_p, sig_fdr, fdr_used

_SUPERSCRIPTS = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")

def format_count(n):
    if n == 0:
        return "0"
    exp = int(math.floor(math.log10(abs(n))))
    if exp < 4:
        return str(n)
    mantissa = n / (10 ** exp)
    exp_str = str(exp).translate(_SUPERSCRIPTS)
    return f"{mantissa:.2g}×10{exp_str}"

def make_plot(stages, counts, pass_labels, out_prefix, bias_overlays=None, trait_order=None, color_map=None, run_label=None):
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    try:
        font_size = 15
        legend_font_size = 14
        span = 0.35
        x = [i * span for i in range(len(stages))]
        heights = [math.log10(c) if c > 0 else math.nan for c in counts]

        fig, ax = plt.subplots(figsize=(8.0, 5.5))
        bar_width = 0.28

        overlay_indices = set()
        if bias_overlays:
            for idx, bdict, tie_ct in bias_overlays:
                total = tie_ct + sum(bdict.values())
                if total > 0:
                    overlay_indices.add(idx)
        base_colors = []
        for i, lbl in enumerate(pass_labels):
            if lbl == "Pass III" or i in overlay_indices:
                base_colors.append("none")
            else:
                base_colors.append("#555555")

        bars = ax.bar(
            x,
            heights,
            color=base_colors,
            edgecolor="#000000",
            linewidth=1.5,
            width=bar_width,
        )
        ax.set_xlim(min(x) - 0.25, max(x) + 0.25)
        ax.set_xticks(x)
        wrapped_stages = []
        for s in stages:
            lbl = s
            if " (" in lbl:
                lbl = lbl.replace(" (", "\n(")
            if lbl.lower().startswith(("unique ", "candidate ", "matching ")):
                lbl = lbl.replace(" ", "\n", 1)
            wrapped_stages.append(lbl)
        ax.set_xticklabels(wrapped_stages, rotation=0, ha="center", fontsize=font_size)
        ax.set_ylabel(r"Feature count (log$_{10}$)", fontsize=font_size)
        title = r"Number of $k$-mers and unitigs after each BiTUGA pass"
        if run_label:
            title += f"\n({run_label})"
        ax.set_title(title, fontsize=font_size + 1, pad=32)
        ax.tick_params(axis="y", labelsize=font_size)

        for i, (rect, count) in enumerate(zip(bars, counts)):
            if count <= 0 or not math.isfinite(rect.get_height()):
                continue
            if pass_labels[i] == "Pass III":
                rect.set_edgecolor("#000000")
                rect.set_linewidth(1.5)
            else:
                rect.set_edgecolor("none")
                rect.set_linewidth(0.0)
            ax.text(
                rect.get_x() + rect.get_width() / 2,
                rect.get_height() + 0.03,
                format_count(count),
                ha="center",
                va="bottom",
                fontsize=font_size,
            )

        legend_handles = []
        if bias_overlays:
            seen_traits = set()
            has_tie = False
            for idx, bias_dict, tie_ct in bias_overlays:
                total = tie_ct + sum(bias_dict.values())
                if total <= 0 or idx >= len(bars):
                    continue
                if pass_labels[idx] == "Pass III":
                    continue
                bar = bars[idx]
                x0 = bar.get_x()
                w = bar.get_width()
                h = bar.get_height()
                y0 = 0.0
                for t in trait_order:
                    val = bias_dict.get(t, 0)
                    if val <= 0:
                        continue
                    ax.add_patch(
                        Rectangle(
                            (x0, y0),
                            w,
                            h * (val / total),
                            color=color_map.get(t, "#555555"),
                            alpha=0.8,
                            linewidth=0,
                        )
                    )
                    y0 += h * (val / total)
                    seen_traits.add(t)
                if tie_ct > 0:
                    ax.add_patch(
                        Rectangle(
                            (x0, y0),
                            w,
                            h * (tie_ct / total),
                            color="#7f8c8d",
                            alpha=0.8,
                            linewidth=0,
                        )
                    )
                    has_tie = True
            
            if seen_traits or has_tie:
                for t in trait_order:
                    if t in seen_traits:
                        legend_handles.append(
                            Rectangle((0, 0), 0.4, 0.4, color=color_map.get(t, "#555555"), alpha=0.8, label=f"{t}-biased")
                        )
                if has_tie:
                    legend_handles.append(Rectangle((0, 0), 0.4, 0.4, color="#7f8c8d", alpha=0.8, label="unbiased"))

        y_min, y_max = ax.get_ylim()
        top = y_max + 2.0
        label_y = 1.005
        for rect, label in zip(bars, pass_labels):
            ax.text(
                rect.get_x() + rect.get_width() / 2,
                label_y,
                label,
                ha="center",
                va="bottom",
                fontsize=font_size,
                color="black",
                transform=ax.get_xaxis_transform(),
                clip_on=False,
            )
        ax.set_ylim(y_min, top)

        if legend_handles:
            ax.legend(handles=legend_handles, fontsize=legend_font_size, frameon=False, loc="upper right", handlelength=0.8, handleheight=0.8)

        fig.tight_layout()
        pdf_path = out_prefix.with_suffix(".pdf")
        png_path = out_prefix.with_suffix(".png")
        sys.stderr.write(f"[DEBUG] Saving plot to {pdf_path} and {png_path}\n")

        fig.savefig(pdf_path, dpi=600)
        fig.savefig(png_path, dpi=600)
        plt.close(fig)
        return True

    except Exception as exc:
        print(f"[ERROR] Plot failed: {exc}", file=sys.stderr)
        traceback.print_exc()
        return False

def write_tsv(stages, counts, out_tsv, bias_info=None):
    def _sanitize(label: str) -> str:
        clean = label.replace("$", "")
        clean = clean.replace(r"\mathrm{raw}", "raw")
        clean = clean.replace(r"\mathrm{adj}", "adj")
        clean = clean.replace(r"\le", "<=")
        clean = clean.replace(r"$k$-mers", "k-mers")
        clean = clean.replace("≤", "<=")
        clean = clean.replace(r"\n", " ")
        clean = clean.replace("\\", "")
        return clean

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    try:
        with out_tsv.open("w") as fh:
            bias_cols = {}
            if bias_info:
                for idx, (bias_dict, tie_ct, traits) in bias_info.items():
                    bias_cols[idx] = (bias_dict, tie_ct, traits)
            fh.write("stage\tcount")
            if bias_cols:
                _, _, traits_example = next(iter(bias_cols.values()))
                for t in traits_example:
                    fh.write(f"\tbias_{t}")
                fh.write("\tbias_ties")
            fh.write("\n")
            for idx, (s, c) in enumerate(zip(stages, counts)):
                fh.write(f"{_sanitize(s)}\t{c}")
                if idx in bias_cols:
                    bdict, ties, traits = bias_cols[idx]
                    for t in traits:
                        fh.write(f"\t{bdict.get(t, 0)}")
                    fh.write(f"\t{ties}")
                fh.write("\n")
        return True
    except Exception as exc:
        print(f"[WARN] Could not write TSV: {exc}", file=sys.stderr)
        return False

def main():
    ap = argparse.ArgumentParser(description="Build feature funnel plot.")
    ap.add_argument("--run-root", required=True, help="Run directory.")
    ap.add_argument("--out-prefix", help="Output prefix.")
    ap.add_argument("--alpha", type=float, default=0.05, help="FDR alpha.")
    ap.add_argument("--p-threshold", type=float, default=0.05, help="Raw p threshold.")
    ap.add_argument(
        "--palette",
        choices=["default", "mint"],
        default="default",
        help="Color palette for traits (default: funnel colors; mint = soft blue/green tones).",
    )
    args = ap.parse_args()

    root = Path(args.run_root).resolve()
    res = root / "results"
    stats_dir = res / "statistics"
    k_val = read_k(root)
    k_label = f"{k_val}-mer" if k_val else r"$k$-mer"

    outdir_default = stats_dir / "plots" / "pipeline_funnel"
    prefix = (
        Path(args.out_prefix).with_suffix("")
        if args.out_prefix
        else outdir_default
    )
    out_plot = prefix
    out_tsv = prefix.with_suffix(".tsv")

    plt, Rectangle = load_matplotlib(out_plot.parent)
    if plt is None:
        sys.exit(0)

    stages = []
    counts = []

    bias_pre, bias_pre_t = parse_pre_bias(stats_dir)
    cand_bias, cand_bias_t = parse_candidate_bias(stats_dir)
    bias_match, bias_match_t = parse_unitig_match_bias(stats_dir)
    fisher_mode, sig_total_summary, bias_sig, bias_sig_t = parse_fisher_summary(stats_dir)

    trait_order = []
    def add_traits(d):
        for t in d:
            if t not in {"tie", "ties", "bias_tie", "bias_ties"} and t not in trait_order:
                trait_order.append(t)
    for d in (bias_pre, cand_bias, bias_match, bias_sig):
        add_traits(d.keys())
    trait_order = sorted(trait_order)
    color_map = get_color_map(trait_order, palette=args.palette)

    union_total = parse_kmer_union(stats_dir)
    stages.append(fr"unique {k_label}s")
    counts.append(union_total)

    cand = res / "candidate_kmers.fasta"
    if not cand.is_file():
        cand_gz = cand.with_suffix(cand.suffix + ".gz")
        cand = cand_gz if cand_gz.is_file() else None
    cand_count = count_fasta_records(cand) if cand else 0
    stages.append(fr"candidate {k_label}s")
    counts.append(cand_count)

    perfect = res / "kmers_perfect_contrast.fasta"
    if not perfect.is_file():
        perfect_gz = perfect.with_suffix(perfect.suffix + ".gz")
        perfect = perfect_gz if perfect_gz.is_file() else None
    if perfect:
        stages.append(r"Perfect $k$-mers")
        counts.append(count_fasta_records(perfect))

    u_stats = parse_unitig_stats(stats_dir)
    stages.append("unitigs")
    counts.append(int(u_stats.get("n_unitigs", 0)))

    stages.append("matching unitigs")
    counts.append(count_presence(res / "unitig_presence.tsv"))

    fdr_enabled = fisher_mode in {"bh", "fdr", "bh_feasible", "storey"}
    detailed_path = stats_dir / "unitig_association_details.tsv"
    tested, sig_p, sig_fdr, fdr_used_calc = parse_fisher(
        detailed_path,
        alpha=args.alpha,
        p_threshold=args.p_threshold,
        fdr_enabled=fdr_enabled,
    )
    if sig_total_summary is not None:
        sig_count = sig_total_summary
    else:
        sig_count = sig_fdr if fdr_enabled else sig_p
    fdr_used = fdr_enabled

    bias_overlay = None
    if fdr_used:
        if fisher_mode == "storey":
            stages.append(rf"unitigs ($q$-value ≤ {args.alpha})")
        else:
            stages.append(rf"unitigs ($p_{{\mathrm{{adj}}}} ≤ {args.alpha}$)")
        counts.append(sig_count)
        bias_overlay = (len(stages) - 1, bias_sig, bias_sig_t)
    else:
        stages.append(rf"unitigs ($p_{{\mathrm{{raw}}}} ≤ {args.p_threshold}$)")
        counts.append(sig_count)
        bias_overlay = (len(stages) - 1, bias_sig, bias_sig_t)

    pass_labels = []
    def stage_is(s, keys):
        for k in keys:
            if k.lower() in s.lower():
                return True
        return False

    for s in stages:
        if stage_is(s, ["unique"]):
            pass_labels.append("Pass I")
        elif stage_is(s, ["candidate", "perfect"]):
            pass_labels.append("Pass II")
        elif stage_is(s, ["matched unitigs", "matching unitigs"]):
            pass_labels.append("Pass IV")
        elif stage_is(s, ["adj", "raw", "q-value", "significant unitigs"]):
            pass_labels.append("Pass V")
        elif stage_is(s, ["unitigs"]):
            pass_labels.append("Pass III")
        else:
            pass_labels.append("")

    overlays = []
    if cand_bias or cand_bias_t > 0:
        idx_cand = next((i for i,s in enumerate(stages) if "candidate" in s.lower()), None)
        if idx_cand is not None:
            overlays.append((idx_cand, cand_bias, cand_bias_t))
    
    if bias_pre or bias_pre_t > 0:
        overlays.append((0, bias_pre, bias_pre_t))
    
    if bias_match or bias_match_t > 0:
        idx_um = next(
            (
                i
                for i, s in enumerate(stages)
                if any(k in s.lower() for k in ("matched unitigs", "matching unitigs"))
            ),
            None,
        )
        if idx_um is not None:
            overlays.append((idx_um, bias_match, bias_match_t))
            
    if bias_overlay:
        overlays.append(bias_overlay)

    bias_info = {}
    for idx, bdict, tie_ct in overlays:
        if not bdict and tie_ct <= 0:
            continue
        if idx not in bias_info:
            cleaned = {k: v for k, v in bdict.items() if k not in {"tie", "ties", "bias_tie", "bias_ties"}}
            bias_info[idx] = (cleaned, tie_ct, trait_order)

    write_tsv(stages, counts, out_tsv, bias_info=bias_info if bias_info else None)

    run_label = root.name.replace("_", " ")

    make_plot(
        stages,
        counts,
        pass_labels,
        out_plot,
        bias_overlays=overlays,
        trait_order=trait_order,
        color_map=color_map,
        run_label=run_label,
    )

    def clean_label(lbl):
        return re.sub(r"[\\$]", "", lbl)

    for s, c in zip(stages, counts):
        print(f"[INFO] {clean_label(s)}: {c}")

if __name__ == "__main__":
    main()
