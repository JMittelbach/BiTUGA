#!/usr/bin/env python3
import argparse
import subprocess
import sys
from collections import defaultdict, OrderedDict
from pathlib import Path
from typing import Optional
import gzip
import re
import tempfile

try:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, ConnectionPatch
    from matplotlib import transforms
    from matplotlib.ticker import FuncFormatter
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
except ImportError:
    plt = None


BWA_INDEX_EXTENSIONS = [".amb", ".ann", ".bwt", ".pac", ".sa"]
BWA_MAX_LENGTH = 100


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def read_chr_lengths(fasta_path: Path) -> OrderedDict:
    chr_lengths = OrderedDict()
    with open_text(fasta_path) as fh:
        chr_name = None
        length = 0
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if chr_name is not None:
                    chr_lengths[chr_name] = length
                chr_name = line[1:].strip().split()[0]
                length = 0
            else:
                length += len(line.strip())
        if chr_name is not None:
            chr_lengths[chr_name] = length
    return chr_lengths


def read_paf_intervals_stream(paf_stream, copy_handle=None):
    intervals_by_chr = defaultdict(list)
    for line in paf_stream:
        line = line.strip()
        if copy_handle is not None:
            copy_handle.write(line + "\n")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        qname = fields[0]
        tname = fields[5]
        try:
            tstart = int(fields[7])
            tend = int(fields[8])
        except ValueError:
            continue
        if tend <= tstart:
            continue
        intervals_by_chr[tname].append((tstart, tend, qname))
    for chr_name in intervals_by_chr:
        intervals_by_chr[chr_name].sort()
    return intervals_by_chr


_cigar_full_match_re = re.compile(r"^(\d+)M$")


def parse_zoom_range(spec: str):
    """
    Parse a zoom range specification in Mb, e.g. "6-12" or "6:12".
    Returns (start_mb, end_mb) as floats, or raises ValueError.
    """
    parts = re.split(r"[,:-]", spec)
    parts = [p for p in parts if p.strip() != ""]
    if len(parts) != 2:
        raise ValueError("Zoom range must be START-END in Mb, e.g. 6-12")
    try:
        start = float(parts[0])
        end = float(parts[1])
    except ValueError:
        raise ValueError("Zoom range bounds must be numeric (Mb)")
    if start < 0 or end <= start:
        raise ValueError("Zoom range must have start >= 0 and end > start (Mb)")
    return start, end


def _save_vector_pdf(fig, out_pdf: Path):
    """
    Save a vector-only PDF alongside the rasterized version by temporarily
    disabling rasterization on artists.
    """
    vector_pdf = out_pdf.with_name(out_pdf.stem + ".vector.pdf")
    raster_states = []
    for ax in fig.axes:
        artists = []
        if hasattr(ax, "collections"):
            artists.extend(ax.collections)
        if hasattr(ax, "lines"):
            artists.extend(ax.lines)
        for artist in artists:
            if hasattr(artist, "get_rasterized") and hasattr(artist, "set_rasterized"):
                raster_states.append((artist, artist.get_rasterized()))
                artist.set_rasterized(False)
    fig.savefig(vector_pdf, dpi=300)
    for artist, state in raster_states:
        artist.set_rasterized(state)
    print(f"[INFO] Saved vector-only PDF to {vector_pdf}")


def md_mismatch_count(md: str) -> int:
    mismatches = 0
    i = 0
    n = len(md)
    while i < n:
        if md[i].isdigit():
            while i < n and md[i].isdigit():
                i += 1
        elif md[i] == "^":
            i += 1
            while i < n and md[i].isalpha():
                mismatches += 1
                i += 1
        else:
            mismatches += 1
            i += 1
    return mismatches


def read_sam_intervals_stream(sam_stream):
    intervals_by_chr = defaultdict(list)
    for line in sam_stream:
        if not line:
            continue
        if line[0] == "@":
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue
        qname = fields[0]
        rname = fields[2]
        if rname == "*" or rname == "":
            continue
        try:
            pos = int(fields[3])
        except ValueError:
            continue
        cigar = fields[5]
        if cigar == "*" or cigar == "":
            continue
        m = _cigar_full_match_re.match(cigar)
        if not m:
            continue
        md = None
        for f in fields[11:]:
            if f.startswith("MD:Z:"):
                md = f[5:]
                break
        if md is None:
            continue
        if md_mismatch_count(md) > 0:
            continue
        read_len = int(m.group(1))
        ref_pos = pos - 1
        start = ref_pos
        end = ref_pos + read_len
        if end > start:
            intervals_by_chr[rname].append((start, end, qname))
    for chr_name in intervals_by_chr:
        intervals_by_chr[chr_name].sort()
    return intervals_by_chr


def compute_window_coverage(chr_len: int, intervals, window_size: int):
    positions = []
    coverage = []
    window_ids = []
    if chr_len <= 0:
        return positions, coverage, window_ids
    n_windows = (chr_len + window_size - 1) // window_size
    i = 0
    for w in range(n_windows):
        win_start = w * window_size
        win_end = min(win_start + window_size, chr_len)
        win_len = win_end - win_start
        if win_len <= 0:
            continue
        while i < len(intervals) and intervals[i][1] <= win_start:
            i += 1
        covered = 0
        j = i
        win_uid_set = set()
        while j < len(intervals) and intervals[j][0] < win_end:
            s, e, uid = intervals[j]
            overlap_start = max(win_start, s)
            overlap_end = min(win_end, e)
            if overlap_end > overlap_start:
                covered += overlap_end - overlap_start
                win_uid_set.add(uid)
            j += 1
        cov_pct = 100.0 * covered / win_len
        positions.append(win_start)
        coverage.append(cov_pct)
        window_ids.append(win_uid_set)
    return positions, coverage, window_ids


def compute_all_window_coverage(chr_lengths: OrderedDict, intervals_by_chr, window_size: int):
    chr_windows = {}
    for chr_name, chr_len in chr_lengths.items():
        intervals = intervals_by_chr.get(chr_name, [])
        positions, cov, ids = compute_window_coverage(chr_len, intervals, window_size)
        chr_windows[chr_name] = (positions, cov, ids)
    return chr_windows


def build_genome_coordinates(chr_lengths: OrderedDict):
    chr_offsets = {}
    chrom_regions = []
    offset = 0
    for chr_name, length in chr_lengths.items():
        chr_offsets[chr_name] = offset
        chrom_regions.append((chr_name, offset, offset + length))
        offset += length
    genome_len = offset
    return chr_offsets, genome_len, chrom_regions


def write_coverage_tsv(chr_lengths: OrderedDict,
                       chr_offsets,
                       window_size: int,
                       chr_windows,
                       out_tsv: Path,
                       panel_label: Optional[str] = None):
    with out_tsv.open("w") as out:
        out.write(
            "chromosome\twindow_start\twindow_end\t"
            "global_start\tglobal_end\tcoverage_pct\n"
        )
        for chr_name, chr_len in chr_lengths.items():
            offset = chr_offsets[chr_name]
            positions, coverages, _ids = chr_windows.get(chr_name, ([], [], []))
            for pos, cov in zip(positions, coverages):
                win_start = pos
                win_end = min(pos + window_size, chr_len)
                global_start = offset + win_start
                global_end = offset + win_end
                out.write(
                    f"{chr_name}\t{win_start}\t{win_end}\t"
                    f"{global_start}\t{global_end}\t{cov:.6f}\n"
                )
    prefix = f"Panel {panel_label}: " if panel_label else ""
    print(f"[INFO] {prefix}Saved coverage table to {out_tsv}")


def plot_coverage_panels(panel_results, out_png: Path, out_pdf: Optional[Path] = None, x_scale: float = 1.0):
    if plt is None:
        print("[ERROR] matplotlib is required for plotting but is not installed.", file=sys.stderr)
        sys.exit(1)

    n_panels = len(panel_results)
    fig_width = 14 * x_scale
    fig, axes = plt.subplots(n_panels, 1, figsize=(fig_width, 3.6 * n_panels), sharex=False)
    if n_panels == 1:
        axes = [axes]

    label_x = -0.035

    for idx, (ax, panel) in enumerate(zip(axes, panel_results)):
        chr_lengths = panel["chr_lengths"]
        chr_offsets = panel["chr_offsets"]
        chrom_regions = panel["chrom_regions"]
        chr_windows = panel["chr_windows"]
        genome_len = panel["genome_len"]

        for chr_name, chr_len in chr_lengths.items():
            positions, cov, _ids = chr_windows.get(chr_name, ([], [], []))
            if not positions:
                continue
            offset = chr_offsets[chr_name]
            pos_global = [p + offset for p in positions]
            ax.plot(
                pos_global,
                cov,
                ".",
                markersize=1,
                alpha=0.7,
                color="black",
                rasterized=True,
            )
        for chrom_idx, (chr_name, start, end) in enumerate(chrom_regions):
            color = "#f0f0f0" if chrom_idx % 2 == 0 else "#ffffff"
            ax.axvspan(start, end, facecolor=color, alpha=0.5, zorder=-1)
        xticks = []
        xticklabels = []
        for chr_name, start, end in chrom_regions:
            center = (start + end) / 2
            xticks.append(center)
            xticklabels.append(chr_name)
        if len(xticks) > 15:
            step = max(2, len(xticks) // 15)
            xticks = xticks[::step]
            xticklabels = xticklabels[::step]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, rotation=0, ha="center", fontsize=11)
        ax.set_xlim(0, genome_len)
        ax.set_ylim(0, 80)
        ax.set_yticks([0, 20, 40, 60, 80])
        ax.set_ylabel("Unitig occupancy (%)", fontsize=12)
        if idx == n_panels - 1:
            ax.set_xlabel("Chromosome", fontsize=12)
        else:
            ax.set_xlabel("")
        ax.tick_params(axis="both", labelsize=11)
        ax.text(
            label_x,
            1.18,
            panel["label"],
            transform=ax.transAxes,
            fontweight="bold",
            fontsize=16,
            ha="left",
            va="top",
        )

    fig.subplots_adjust(left=0.06, right=0.99, top=0.95, bottom=0.08, hspace=0.24)
    fig.savefig(out_png, dpi=1200)
    print(f"[INFO] Saved plot to {out_png}")
    if out_pdf is not None:
        fig.savefig(out_pdf, dpi=1200)
        print(f"[INFO] Saved plot to {out_pdf}")
        _save_vector_pdf(fig, out_pdf)


def plot_coverage_with_zooms(panel_results,
                             zoom_targets,
                             zoom_ranges,
                             cut_ranges,
                             out_png: Path,
                             out_pdf: Optional[Path] = None,
                             x_scale: float = 1.0,
                             insert_mode: bool = False,
                             zoom_fill: bool = False):
    if plt is None:
        print("[ERROR] matplotlib is required for plotting but is not installed.", file=sys.stderr)
        sys.exit(1)
    n_rows = len(panel_results)
    fig_width = 14 * x_scale
    if insert_mode:
        fig, axes = plt.subplots(n_rows, 1, figsize=(fig_width, 3.6 * n_rows), sharex=False)
        if n_rows == 1:
            axes = [axes]
    else:
        fig, axes = plt.subplots(
            n_rows,
            2,
            figsize=(fig_width, 3.6 * n_rows),
            sharex=False,
            gridspec_kw={"width_ratios": [4, 1]},
        )
        if n_rows == 1:
            axes = [axes]

    label_main_x = -0.045
    label_zoom_x = -0.05
    label_iter = (chr(ord("A") + i) for i in range(26))

    for idx, panel in enumerate(panel_results):
        row_axes = axes[idx] if insert_mode else axes[idx]
        if insert_mode:
            ax_main = row_axes
            ax_zoom = None
            main_label = chr(ord("A") + idx)
        else:
            ax_main, ax_zoom = row_axes
            main_label = next(label_iter, "A")
        chr_lengths = panel["chr_lengths"]
        chr_offsets = panel["chr_offsets"]
        chrom_regions = panel["chrom_regions"]
        chr_windows = panel["chr_windows"]
        genome_len = panel["genome_len"]

        for chr_name, chr_len in chr_lengths.items():
            positions, cov, _ids = chr_windows.get(chr_name, ([], [], []))
            if not positions:
                continue
            offset = chr_offsets[chr_name]
            pos_global = [p + offset for p in positions]
            ax_main.plot(pos_global, cov, ".", markersize=1, alpha=1.0, color="black", rasterized=True)
        for chrom_idx, (chr_name, start, end) in enumerate(chrom_regions):
            color = "#f0f0f0" if chrom_idx % 2 == 0 else "#ffffff"
            ax_main.axvspan(start, end, facecolor=color, alpha=0.5, zorder=-1)
        xticks = []
        xticklabels = []
        for chr_name, start, end in chrom_regions:
            center = (start + end) / 2
            xticks.append(center)
            xticklabels.append(chr_name)
        if len(xticks) > 15:
            step = max(2, len(xticks) // 15)
            xticks = xticks[::step]
            xticklabels = xticklabels[::step]
        ax_main.set_xticks(xticks)
        ax_main.set_xticklabels(xticklabels, rotation=0, ha="center", fontsize=11)
        ax_main.set_xlim(0, genome_len)
        ax_main.set_ylim(0, 80)
        ax_main.set_yticks([0, 20, 40, 60, 80])
        ax_main.set_ylabel("Unitig occupancy (%)", fontsize=12)
        if idx == n_rows - 1:
            ax_main.set_xlabel("Chromosome", fontsize=12)
        else:
            ax_main.set_xlabel("")
        ax_main.tick_params(axis="both", labelsize=11)
        ax_main.text(
            label_main_x,
            1.10,
            main_label,
            transform=ax_main.transAxes,
            fontweight="bold",
            fontsize=16,
            ha="left",
            va="top",
        )
        if panel["label"] == "A":
            ax_main.set_title("P. tremula", fontstyle="italic", fontsize=13, loc="left", pad=10)
        elif panel["label"] == "B":
            ax_main.set_title("G. biloba", fontstyle="italic", fontsize=13, loc="left", pad=10)

        zoom_chr = zoom_targets.get(panel["label"])
        zoom_range_mb = zoom_ranges.get(panel["label"])
        cut_range_mb = cut_ranges.get(panel["label"]) if cut_ranges else None
        if zoom_chr:
            out_zoom_png, out_zoom_pdf = make_zoom_paths(panel["prefix"], zoom_chr)
            if insert_mode:
                if zoom_chr in chr_offsets:
                    z_start = chr_offsets[zoom_chr]
                    z_end = chr_offsets[zoom_chr] + chr_lengths[zoom_chr]
                    if zoom_range_mb:
                        z_start = z_start + zoom_range_mb[0] * 1_000_000
                        z_end = chr_offsets[zoom_chr] + min(chr_lengths[zoom_chr], zoom_range_mb[1] * 1_000_000)
                    y0, y1 = ax_main.get_ylim()
                    pad = (y1 - y0) * 0.02
                    y_start = y0 + pad
                    y_height = max((y1 - y0) - 2 * pad, (y1 - y0) * 0.5)
                    if zoom_fill:
                        ax_main.axvspan(
                            z_start,
                            z_end,
                            facecolor="#c8f5df",
                            alpha=0.35,
                            zorder=0.5,
                        )
                    else:
                        rect = Rectangle(
                            (z_start, y_start),
                            z_end - z_start,
                            y_height,
                            linewidth=1,
                            edgecolor="gray" if not zoom_fill else "black",
                            facecolor="none",
                            linestyle="--",
                            zorder=5,
                        )
                        ax_main.add_patch(rect)
                w_inset = 0.24
                h_inset = 0.45
                x0, x1 = ax_main.get_xlim()
                start_norm = (z_start - x0) / (x1 - x0) if x1 > x0 else 0.5
                end_norm = (z_end - x0) / (x1 - x0) if x1 > x0 else 0.5
                center_norm = (start_norm + end_norm) / 2
                if panel["label"] == "B":
                    x_box = min(max(end_norm + 0.05, 0.08), 1 - w_inset - 0.05)
                else:
                    x_box = min(max(center_norm - w_inset / 2, 0.05), 1 - w_inset - 0.05)
                y_box = 0.90 - h_inset
                inset = inset_axes(
                    ax_main,
                    width="100%",
                    height="100%",
                    loc="lower left",
                    bbox_to_anchor=(x_box, y_box, w_inset, h_inset),
                    bbox_transform=ax_main.transAxes,
                    borderpad=0.1,
                )
                plot_zoom_chromosome(
                    panel,
                    zoom_chr,
                    out_zoom_png,
                    out_zoom_pdf,
                    zoom_range_mb=zoom_range_mb,
                    ax=inset,
                    show_ylabel=False,
                    cut_range_mb=cut_range_mb,
                )
                try:
                    rect_center = (
                        z_end if panel["label"] == "B" else z_start,
                        y_start + y_height / 2,
                    )
                    inset_anchor = (0.0, 0.0) if panel["label"] == "B" else (1.0, 0.0)
                    conn = ConnectionPatch(
                        xyA=rect_center,
                        coordsA=ax_main.transData,
                        xyB=inset_anchor,
                        coordsB=inset.transAxes,
                        linestyle="--",
                        color="black",
                        linewidth=0.8,
                        zorder=1,
                    )
                    ax_main.add_artist(conn)
                except Exception:
                    pass
            else:
                zoom_label = next(label_iter, "A")
                plot_zoom_chromosome(
                    panel,
                    zoom_chr,
                    out_zoom_png,
                    out_zoom_pdf,
                    zoom_range_mb=zoom_range_mb,
                    ax=ax_zoom,
                    show_ylabel=False,
                    cut_range_mb=cut_range_mb,
                )
                ax_zoom.text(
                    label_zoom_x,
                    1.08,
                    zoom_label,
                    transform=ax_zoom.transAxes,
                    fontweight="bold",
                    fontsize=16,
                    ha="left",
                    va="top",
                )
        elif not insert_mode:
            ax_zoom.axis("off")

    fig.subplots_adjust(left=0.06, right=0.985, top=0.95, bottom=0.08, hspace=0.24, wspace=0.08)
    fig.savefig(out_png, dpi=1200)
    print(f"[INFO] Saved plot with zooms to {out_png}")
    if out_pdf is not None:
        fig.savefig(out_pdf, dpi=1200)
        print(f"[INFO] Saved plot with zooms to {out_pdf}")
        _save_vector_pdf(fig, out_pdf)


def run_minimap2(minimap2_bin: str, ref: Path, unitigs: Path, threads: int):
    cmd = [minimap2_bin, "-t", str(threads), str(ref), str(unitigs)]
    print(f"[INFO] Running: {' '.join(cmd)}", file=sys.stderr)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    return proc

def run_minimap2_strict(minimap2_bin: str, ref: Path, unitigs: Path, threads: int):
    cmd = [
        minimap2_bin,
        "-x", "asm10",       
        "-t", str(threads),
        str(ref),
        str(unitigs),
    ]
    print(f"[INFO] Running (asm10 mode): {' '.join(cmd)}", file=sys.stderr)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    return proc


def _bwa_index_complete(prefix: Path) -> bool:
    return all(Path(f"{prefix}{ext}").exists() for ext in BWA_INDEX_EXTENSIONS)


def find_existing_bwa_index(ref: Path) -> Optional[Path]:
    """
    Check for an existing BWA index alongside the reference. We accept both
    .fa.* and .fasta.* index prefixes in the same directory.
    """
    candidates = [ref]
    suffix = ref.suffix.lower()
    if suffix in (".fa", ".fasta"):
        alt_suffix = ".fasta" if suffix == ".fa" else ".fa"
        candidates.append(ref.with_suffix(alt_suffix))
    for candidate in candidates:
        if _bwa_index_complete(candidate):
            return candidate
    return None


def ensure_bwa_index(bwa_bin: str, ref: Path) -> Path:
    existing_prefix = find_existing_bwa_index(ref)
    if existing_prefix is not None:
        print(f"[INFO] BWA index found for {existing_prefix}, skipping indexing", file=sys.stderr)
        return existing_prefix

    print(f"[INFO] BWA index not found for {ref}, running 'bwa index'...", file=sys.stderr)
    cmd = [bwa_bin, "index", str(ref)]
    ret = subprocess.run(cmd).returncode
    if ret != 0:
        print(f"[ERROR] bwa index exited with code {ret}", file=sys.stderr)
        sys.exit(ret)
    print("[INFO] BWA index created", file=sys.stderr)
    return ref


def run_bwa_aln_samse(bwa_bin: str, ref: Path, unitigs: Path, threads: int, seed_length: int):
    index_prefix = ensure_bwa_index(bwa_bin, ref)
    tmpdir = tempfile.TemporaryDirectory()
    sai_path = Path(tmpdir.name) / "unitigs.sai"
    aln_cmd = [
        bwa_bin,
        "aln",
        "-t", str(threads),
        "-n", "0",
        "-o", "0",
        "-e", "0",
        "-l", str(seed_length),
        str(index_prefix),
        str(unitigs),
    ]
    print(f"[INFO] Running: {' '.join(aln_cmd)}", file=sys.stderr)
    with sai_path.open("w") as sai_out:
        ret = subprocess.run(aln_cmd, stdout=sai_out).returncode
    if ret != 0:
        print(f"[ERROR] bwa aln exited with code {ret}", file=sys.stderr)
        tmpdir.cleanup()
        sys.exit(ret)
    samse_cmd = [bwa_bin, "samse", str(index_prefix), str(sai_path), str(unitigs)]
    print(f"[INFO] Running: {' '.join(samse_cmd)}", file=sys.stderr)
    proc = subprocess.Popen(samse_cmd, stdout=subprocess.PIPE, text=True)
    proc._bwa_tmpdir = tmpdir
    return proc


def make_panel_prefix(unitigs: Path, mode: str, window: int, outdir: Path, label: str, use_label_prefix: bool):
    base = f"{unitigs.name}.{mode}.win{window}"
    if use_label_prefix:
        base = f"{label}.{base}"
    return outdir / base


def process_panel(label: str,
                  ref: Path,
                  unitigs: Path,
                  outdir: Path,
                  window: int,
                  minimap2_bin: str,
                  bwa_bin: str,
                  threads: int,
                  mode: str,
                  top_windows: int,
                  export_chr,
                  bwa_seed_length: int,
                  male_biased_only: bool,
                  use_label_prefix: bool):
    if not ref.exists():
        print(f"[ERROR] Panel {label}: Reference file does not exist: {ref}", file=sys.stderr)
        sys.exit(1)
    if not unitigs.exists():
        print(f"[ERROR] Panel {label}: Unitig file does not exist: {unitigs}", file=sys.stderr)
        sys.exit(1)
    if str(ref).endswith(".gz") and mode in ("bwa", "hybrid"):
        print(f"[ERROR] Panel {label}: BWA and hybrid modes require uncompressed reference FASTA", file=sys.stderr)
        sys.exit(1)

    tmpdir_male = None
    unitigs_to_use = unitigs
    if male_biased_only:
        print(f"[INFO] Panel {label}: Filtering to male-biased unitigs (male > female)", file=sys.stderr)
        tmpdir_male = tempfile.TemporaryDirectory()
        unitigs_to_use = filter_male_biased(unitigs, tmpdir_male)
        if count_unitigs(unitigs_to_use) == 0:
            print(f"[WARN] Panel {label}: No male-biased unitigs found; outputs will be empty.", file=sys.stderr)

    prefix = make_panel_prefix(unitigs_to_use, mode, window, outdir, label, use_label_prefix)
    out_tsv = prefix.with_suffix(".tsv")

    print(f"[INFO] Panel {label}: Mode {mode}", file=sys.stderr)
    if mode == "hybrid":
        print(
            f"[INFO] Panel {label}: BWA for unitigs <= {BWA_MAX_LENGTH} bp "
            f"(exact matches, seed length {bwa_seed_length})",
            file=sys.stderr,
        )
        print(f"[INFO] Panel {label}: minimap2 for unitigs > {BWA_MAX_LENGTH} bp", file=sys.stderr)
    elif mode == "bwa":
        print(
            f"[INFO] Panel {label}: BWA for all unitigs "
            f"(exact matches, seed length {bwa_seed_length})",
            file=sys.stderr,
        )
    else:
        print(f"[INFO] Panel {label}: minimap2 for all unitigs", file=sys.stderr)
    print(f"[INFO] Panel {label}: Output TSV: {out_tsv}", file=sys.stderr)

    chr_lengths = read_chr_lengths(ref)
    if not chr_lengths:
        print(f"[ERROR] Panel {label}: Failed to read any chromosomes from reference.", file=sys.stderr)
        sys.exit(1)

    intervals_by_chr = defaultdict(list)

    if mode == "minimap":
        paf_path = prefix.with_suffix(".paf")
        proc = run_minimap2(minimap2_bin, ref, unitigs_to_use, threads)
        try:
            with paf_path.open("w") as paf_out:
                intervals = read_paf_intervals_stream(proc.stdout, copy_handle=paf_out)
        finally:
            ret = proc.wait()
        if ret != 0:
            print(f"[ERROR] Panel {label}: minimap2 exited with code {ret}", file=sys.stderr)
            sys.exit(ret)
        print(f"[INFO] Panel {label}: Saved PAF to {paf_path}", file=sys.stderr)
        intervals_by_chr = merge_intervals(intervals_by_chr, intervals)
    elif mode == "minimap-strict":
        paf_path = prefix.with_suffix(".paf")
        proc = run_minimap2_strict(minimap2_bin, ref, unitigs_to_use, threads)
        try:
            with paf_path.open("w") as paf_out:
                intervals = read_paf_intervals_stream(proc.stdout, copy_handle=paf_out)
        finally:
            ret = proc.wait()
        if ret != 0:
            print(f"[ERROR] Panel {label}: minimap2 exited with code {ret}", file=sys.stderr)
            sys.exit(ret)
        print(f"[INFO] Panel {label}: Saved PAF to {paf_path}", file=sys.stderr)
        intervals_by_chr = merge_intervals(intervals_by_chr, intervals)
    elif mode == "bwa":
        proc = run_bwa_aln_samse(bwa_bin, ref, unitigs_to_use, threads, bwa_seed_length)
        try:
            intervals = read_sam_intervals_stream(proc.stdout)
        finally:
            ret = proc.wait()
            if hasattr(proc, "_bwa_tmpdir"):
                proc._bwa_tmpdir.cleanup()
        if ret != 0:
            print(f"[ERROR] Panel {label}: bwa exited with code {ret}", file=sys.stderr)
            sys.exit(ret)
        intervals_by_chr = merge_intervals(intervals_by_chr, intervals)
    else:
        tmpdir_unitigs = tempfile.TemporaryDirectory()
        try:
            short_fa, long_fa, has_short, has_long = split_unitigs_by_length(
                unitigs_to_use, BWA_MAX_LENGTH, Path(tmpdir_unitigs.name)
            )
            if has_long:
                paf_path = prefix.with_suffix(".paf")
                proc_long = run_minimap2(minimap2_bin, ref, long_fa, threads)
                try:
                    with paf_path.open("w") as paf_out:
                        intervals_long = read_paf_intervals_stream(proc_long.stdout, copy_handle=paf_out)
                finally:
                    ret_long = proc_long.wait()
                if ret_long != 0:
                    print(f"[ERROR] Panel {label}: minimap2 exited with code {ret_long}", file=sys.stderr)
                    sys.exit(ret_long)
                print(f"[INFO] Panel {label}: Saved PAF to {paf_path}", file=sys.stderr)
                intervals_by_chr = merge_intervals(intervals_by_chr, intervals_long)
            if has_short:
                proc_short = run_bwa_aln_samse(
                    bwa_bin, ref, short_fa, threads, bwa_seed_length
                )
                try:
                    intervals_short = read_sam_intervals_stream(proc_short.stdout)
                finally:
                    ret_short = proc_short.wait()
                    if hasattr(proc_short, "_bwa_tmpdir"):
                        proc_short._bwa_tmpdir.cleanup()
                if ret_short != 0:
                    print(f"[ERROR] Panel {label}: bwa exited with code {ret_short}", file=sys.stderr)
                    sys.exit(ret_short)
                intervals_by_chr = merge_intervals(intervals_by_chr, intervals_short)
        finally:
            tmpdir_unitigs.cleanup()

    if not intervals_by_chr:
        print(f"[WARN] Panel {label}: No intervals read from mapper output; plot and table will be empty.", file=sys.stderr)

    if export_chr:
        chr_unitigs = collect_unitigs_for_chromosomes(intervals_by_chr, export_chr)
        export_unitigs_to_fasta(unitigs_to_use, chr_unitigs, prefix)

    chr_windows = compute_all_window_coverage(chr_lengths, intervals_by_chr, window)
    chr_offsets, genome_len, chrom_regions = build_genome_coordinates(chr_lengths)

    write_coverage_tsv(chr_lengths, chr_offsets, window, chr_windows, out_tsv, panel_label=label)

    total_unitigs = count_unitigs(unitigs_to_use)
    lengths_map, headers_map, total_unitig_bases = read_unitig_metadata(unitigs_to_use)
    mapping_counts = defaultdict(int)
    for ints in intervals_by_chr.values():
        for _s, _e, uid in ints:
            mapping_counts[uid] += 1
    mapped_ids = set(mapping_counts.keys())
    total_mapped_unitigs = len(mapped_ids)
    total_mapped_bases = sum(lengths_map.get(uid, 0) for uid in mapped_ids)
    multimapper_ids = {uid for uid, cnt in mapping_counts.items() if cnt > 1}

    if top_windows > 0:
        top_n = top_windows
        top_entries = []
        for chr_name, chr_len in chr_lengths.items():
            positions, coverages, ids_per_window = chr_windows.get(chr_name, ([], [], []))
            offset = chr_offsets[chr_name]
            for pos, cov, win_ids in zip(positions, coverages, ids_per_window):
                win_start = pos
                win_end = min(pos + window, chr_len)
                g_start = offset + win_start
                g_end = offset + win_end
                top_entries.append((cov, chr_name, win_start, win_end, g_start, g_end, win_ids))
        top_entries.sort(key=lambda t: t[0], reverse=True)
        top_entries_all = top_entries
        top_entries = top_entries_all[:top_n]
        top_entries_50 = top_entries_all[:50]
        out_top = prefix.with_suffix(".top_windows.tsv")
        with out_top.open("w") as fh:
            fh.write("#chromosome_stats\tcount_mapped\tpct_of_unitigs\tnorm_unitigs_per_mb\n")
            for chr_name, chr_len in chr_lengths.items():
                ids_chr = {uid for _s, _e, uid in intervals_by_chr.get(chr_name, [])}
                count_chr = len(ids_chr)
                pct = 100.0 * count_chr / total_unitigs if total_unitigs > 0 else 0.0
                norm = count_chr / (chr_len / 1_000_000) if chr_len > 0 else 0.0
                fh.write(f"#chr_stats\t{chr_name}\t{count_chr}\t{pct:.6f}\t{norm:.6f}\n")
            fh.write(f"#total_unitigs\t{total_unitigs}\n")
            fh.write(f"#mapping_rate_unitigs_pct\t{100.0 * total_mapped_unitigs / total_unitigs if total_unitigs > 0 else 0.0:.6f}\n")
            fh.write(f"#multimapper_unitigs\t{len(multimapper_ids)}\n")
            fh.write(f"#multimapper_pct_total\t{100.0 * len(multimapper_ids) / total_unitigs if total_unitigs > 0 else 0.0:.6f}\n")
            fh.write(f"#multimapper_pct_mapped\t{100.0 * len(multimapper_ids) / total_mapped_unitigs if total_mapped_unitigs > 0 else 0.0:.6f}\n")
            fh.write("#chromosome_stats_mapped_only\tcount_mapped\tpct_of_mapped_unitigs\tnorm_mapped_unitigs_per_mb\tlen_norm_pct_of_mapped\n")
            mapped_rates = {}
            total_rate = 0.0
            for chr_name, chr_len in chr_lengths.items():
                ids_chr = {uid for _s, _e, uid in intervals_by_chr.get(chr_name, [])}
                rate = count_chr = len(ids_chr)
                rate = count_chr / (chr_len / 1_000_000) if chr_len > 0 else 0.0
                mapped_rates[chr_name] = rate
                total_rate += rate
            for chr_name, chr_len in chr_lengths.items():
                ids_chr = {uid for _s, _e, uid in intervals_by_chr.get(chr_name, [])}
                count_chr = len(ids_chr)
                pct_mapped = 100.0 * count_chr / total_mapped_unitigs if total_mapped_unitigs > 0 else 0.0
                norm_mapped = count_chr / (chr_len / 1_000_000) if chr_len > 0 else 0.0
                len_norm_pct = 100.0 * mapped_rates.get(chr_name, 0.0) / total_rate if total_rate > 0 else 0.0
                fh.write(f"#chr_stats_mapped\t{chr_name}\t{count_chr}\t{pct_mapped:.6f}\t{norm_mapped:.6f}\t{len_norm_pct:.6f}\n")
            fh.write(f"#total_mapped_unitigs\t{total_mapped_unitigs}\n")
            fh.write("#chromosome_base_stats\tbases_mapped\tpct_of_bases\tnorm_bases_per_mb\n")
            for chr_name, chr_len in chr_lengths.items():
                ids_chr = {uid for _s, _e, uid in intervals_by_chr.get(chr_name, [])}
                bases_chr = sum(lengths_map.get(uid, 0) for uid in ids_chr)
                pct_bases = 100.0 * bases_chr / total_unitig_bases if total_unitig_bases > 0 else 0.0
                norm_bases = bases_chr / (chr_len / 1_000_000) if chr_len > 0 else 0.0
                fh.write(f"#chr_base_stats\t{chr_name}\t{bases_chr}\t{pct_bases:.6f}\t{norm_bases:.6f}\n")
            fh.write(f"#total_bases\t{total_unitig_bases}\n")
            fh.write(f"#mapping_rate_bases_pct\t{100.0 * total_mapped_bases / total_unitig_bases if total_unitig_bases > 0 else 0.0:.6f}\n")
            fh.write("#chromosome_base_stats_mapped_only\tbases_mapped\tpct_of_mapped_bases\n")
            for chr_name, chr_len in chr_lengths.items():
                ids_chr = {uid for _s, _e, uid in intervals_by_chr.get(chr_name, [])}
                bases_chr = sum(lengths_map.get(uid, 0) for uid in ids_chr)
                pct_bases_mapped = 100.0 * bases_chr / total_mapped_bases if total_mapped_bases > 0 else 0.0
                fh.write(f"#chr_base_stats_mapped\t{chr_name}\t{bases_chr}\t{pct_bases_mapped:.6f}\n")
            fh.write(f"#total_mapped_bases\t{total_mapped_bases}\n")
            chr_targets = {"A": "19", "B": "2"}
            target_chr = chr_targets.get(label)
            if target_chr:
                longest = longest_unitig_for_chr(intervals_by_chr, target_chr, lengths_map, headers_map, require_male_biased=True)
                if longest:
                    uid, length_bp, start, end, pval, header_text = longest
                    fh.write(f"#longest_chr{target_chr}\t{uid}\t{length_bp}\t{start}\t{end}\t{pval}\t{header_text}\n")
            def summarize_top(entries):
                if not entries:
                    return None
                chr_counter = defaultdict(int)
                chr_minmax = {}
                for _cov, chr_name, w_start, w_end, _gs, _ge, _ids in entries:
                    chr_counter[chr_name] += 1
                    if chr_name not in chr_minmax:
                        chr_minmax[chr_name] = [w_start, w_end]
                    else:
                        chr_minmax[chr_name][0] = min(chr_minmax[chr_name][0], w_start)
                        chr_minmax[chr_name][1] = max(chr_minmax[chr_name][1], w_end)
                top_chr = max(chr_counter.items(), key=lambda kv: kv[1])[0]
                pct_top = 100.0 * chr_counter[top_chr] / len(entries)
                min_start, max_end = chr_minmax[top_chr]
                return top_chr, pct_top, min_start, max_end

            info_top = summarize_top(top_entries)
            if info_top:
                top_chr, pct_top, min_start, max_end = info_top
                fh.write(f"#top_chr_top_windows\t{top_chr}\t{pct_top:.2f}\t{min_start}\t{max_end}\n")
            info_top50 = summarize_top(top_entries_50)
            if info_top50:
                top_chr, pct_top, min_start, max_end = info_top50
                fh.write(f"#top_chr_top50_windows\t{top_chr}\t{pct_top:.2f}\t{min_start}\t{max_end}\n")
            fh.write("chromosome\twindow_start\twindow_end\tglobal_start\tglobal_end\tcoverage_pct\tunitigs\n")
            for cov, chr_name, win_start, win_end, g_start, g_end, win_ids in top_entries:
                fh.write(f"{chr_name}\t{win_start}\t{win_end}\t{g_start}\t{g_end}\t{cov:.6f}\t{','.join(sorted(win_ids))}\n")
        print(f"[INFO] Panel {label}: Saved top {len(top_entries)} windows to {out_top}", file=sys.stderr)
        top_ids = set()
        for _cov, _chr_name, _wstart, _wend, _gs, _ge, win_ids in top_entries:
            top_ids.update(win_ids)
        if top_ids:
            out_top_fa = prefix.with_suffix(".top_windows.unitigs.fa")
            write_unitigs_subset(unitigs_to_use, top_ids, out_top_fa)

    chr_targets = {"A": "19", "B": "2"}
    target_chr = chr_targets.get(label)
    if target_chr:
        ids_chr = {uid for _s, _e, uid in intervals_by_chr.get(target_chr, [])}
        if ids_chr:
            out_chr_fa = prefix.with_suffix(f".chr{target_chr}.unitigs.fa")
            write_unitigs_subset(unitigs_to_use, ids_chr, out_chr_fa)

    if tmpdir_male:
        tmpdir_male.cleanup()

    return {
        "label": label,
        "ref": ref,
        "unitigs": unitigs,
        "prefix": prefix,
        "chr_lengths": chr_lengths,
        "chr_offsets": chr_offsets,
        "chrom_regions": chrom_regions,
        "chr_windows": chr_windows,
        "window": window,
        "genome_len": genome_len,
    }


def make_zoom_paths(prefix: Path, chrom: str):
    safe_chr = chrom.replace("/", "_")
    base = f"{prefix}.chr{safe_chr}.zoom"
    return Path(f"{base}.png"), Path(f"{base}.pdf")


def plot_zoom_chromosome(panel,
                         chrom_name: str,
                         out_png: Optional[Path] = None,
                         out_pdf: Optional[Path] = None,
                         zoom_range_mb=None,
                         ax=None,
                         show_ylabel: bool = True,
                         cut_range_mb=None):
    if plt is None:
        print("[ERROR] matplotlib is required for plotting but is not installed.", file=sys.stderr)
        sys.exit(1)
    
    chr_windows = panel["chr_windows"]
    chr_lengths = panel["chr_lengths"]
    
    if chrom_name not in chr_lengths:
        print(f"[WARN] Panel {panel['label']}: Chromosome {chrom_name} not found; zoom plot skipped.", file=sys.stderr)
        return
        
    positions, cov, _ids = chr_windows.get(chrom_name, ([], [], []))
    if not positions:
        print(f"[WARN] Panel {panel['label']}: No coverage data for chromosome {chrom_name}; zoom plot skipped.", file=sys.stderr)
        return
        
    xs_mb = [p / 1_000_000 for p in positions]
    chr_len_mb = chr_lengths[chrom_name] / 1_000_000
    
    if zoom_range_mb:
        z_start, z_end = zoom_range_mb
    else:
        z_start, z_end = (0, chr_len_mb)

    created_fig = ax is None
    if created_fig:
        fig, ax = plt.subplots(figsize=(10, 4))
    else:
        fig = ax.figure

    lw = 1.0 


    if cut_range_mb:
        c_start, c_end = cut_range_mb
        
        if c_start > z_start and c_end < z_end:
            ax.axis('off')
            
            ax.plot([0, 1], [1, 1], transform=ax.transAxes, 
                    color='black', linewidth=lw, clip_on=False)

            len_left = c_start - z_start
            len_right = z_end - c_end
            total_vis = len_left + len_right
            gap = 0.02 
            width_avail = 1.0 - gap
            
            w_left = (len_left / total_vis) * width_avail
            w_right = (len_right / total_vis) * width_avail
            
            ax1 = ax.inset_axes([0, 0, w_left, 1.0])
            ax2 = ax.inset_axes([w_left + gap, 0, w_right, 1.0])
            
            idxs_l = [i for i, x in enumerate(xs_mb) if z_start <= x <= c_start]
            ax1.plot([xs_mb[i] for i in idxs_l], [cov[i] for i in idxs_l], 
                     ".", markersize=1, alpha=1.0, color="black", rasterized=True)
            
            ax1.set_xlim(z_start, c_start)
            ax1.set_ylim(0, 80)
            
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.spines['left'].set_visible(True)
            ax1.spines['bottom'].set_visible(True)
            
            ax1.spines['left'].set_linewidth(lw)
            ax1.spines['bottom'].set_linewidth(lw)
            
            ax1.yaxis.tick_left()
            if show_ylabel:
                ax1.set_ylabel(f"Chr {chrom_name} Occupancy (%)", fontsize=12)

            if panel.get("label") == "B":
                ax1.set_xticks([190, 200, 210])

            idxs_r = [i for i, x in enumerate(xs_mb) if c_end <= x <= z_end]
            ax2.plot([xs_mb[i] for i in idxs_r], [cov[i] for i in idxs_r], 
                     ".", markersize=1, alpha=1.0, color="black", rasterized=True)
            
            ax2.set_xlim(c_end, z_end)
            ax2.set_ylim(0, 80)
            
            ax2.spines['top'].set_visible(False)
            ax2.spines['left'].set_visible(False)
            ax2.spines['right'].set_visible(True)
            ax2.spines['bottom'].set_visible(True)
            
            ax2.spines['right'].set_linewidth(lw)
            ax2.spines['bottom'].set_linewidth(lw)

            ax2.yaxis.tick_right()
            ax2.set_yticks([]) 
            
            if panel.get("label") == "B":
                ax2.set_xticks([280, 290, 300])

            d = 0.015 
            kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False, linewidth=lw)
            ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs) 

            kwargs.update(transform=ax2.transAxes)
            ax2.plot((-d, +d), (-d, +d), **kwargs)

            ax.text(0.5, -0.25, "Mb", transform=ax.transAxes, 
                    ha='center', va='top', fontsize=12)
            
            for sub_ax in (ax1, ax2):
                sub_ax.tick_params(axis="both", labelsize=11, width=lw)
                sub_ax.set_yticks([0, 40, 80])

            if created_fig:
                if out_png: fig.savefig(out_png, dpi=1200, bbox_inches='tight')
                if out_pdf: fig.savefig(out_pdf, dpi=1200, bbox_inches='tight')
            
            return

    idxs = [i for i, x in enumerate(xs_mb) if z_start <= x <= z_end]
    ax.plot([xs_mb[i] for i in idxs], [cov[i] for i in idxs], 
            ".", markersize=1, alpha=1.0, color="black", rasterized=True)
    
    ax.set_xlim(z_start, z_end)
    ax.set_ylim(0, 80)
    
    for sp in ax.spines.values():
        sp.set_linewidth(lw)
    
    ax.set_xlabel("Mb", fontsize=12)
    if show_ylabel:
        ax.set_ylabel(f"Chr {chrom_name} Unitig occupancy (%)", fontsize=12)
    
    ax.tick_params(axis="both", labelsize=11, width=lw)
    ax.set_yticks([0, 40, 80])

    if created_fig:
        fig.tight_layout()
        if out_png:
            fig.savefig(out_png, dpi=1200)
            print(f"[INFO] Saved zoom plot for {chrom_name} to {out_png}", file=sys.stderr)
        if out_pdf:
            fig.savefig(out_pdf, dpi=1200)

def split_unitigs_by_length(unitigs: Path, threshold: int, tmpdir: Path):
    short_path = tmpdir / "unitigs_short.fa"
    long_path = tmpdir / "unitigs_long.fa"
    wrote_short = False
    wrote_long = False
    n_short = 0
    n_long = 0
    with open_text(unitigs) as fin, short_path.open("w") as fout_short, long_path.open("w") as fout_long:
        header = None
        seq_lines = []
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_lines).strip()
                    if seq:
                        if len(seq) <= threshold:
                            fout_short.write(f">{header}\n{seq}\n")
                            wrote_short = True
                            n_short += 1
                        else:
                            fout_long.write(f">{header}\n{seq}\n")
                            wrote_long = True
                            n_long += 1
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            seq = "".join(seq_lines).strip()
            if seq:
                if len(seq) <= threshold:
                    fout_short.write(f">{header}\n{seq}\n")
                    wrote_short = True
                    n_short += 1
                else:
                    fout_long.write(f">{header}\n{seq}\n")
                    wrote_long = True
                    n_long += 1
    print(f"[INFO] Unitigs <= {threshold} bp (BWA): {n_short}", file=sys.stderr)
    print(f"[INFO] Unitigs >  {threshold} bp (minimap2): {n_long}", file=sys.stderr)
    return short_path, long_path, wrote_short, wrote_long


def merge_intervals(dict_a, dict_b):
    if not dict_a:
        merged = defaultdict(list)
    else:
        merged = defaultdict(list, {k: v[:] for k, v in dict_a.items()})
    for chr_name, ints in dict_b.items():
        merged[chr_name].extend(ints)
    for chr_name in merged:
        merged[chr_name].sort()
    return merged


def collect_unitigs_for_chromosomes(intervals_by_chr, chromosomes):
    target = OrderedDict()
    for chr_name in chromosomes:
        if chr_name in target:
            continue
        ids = {uid for _start, _end, uid in intervals_by_chr.get(chr_name, [])}
        if not ids:
            print(f"[WARN] No unitigs mapped to {chr_name}; skipping FASTA export.", file=sys.stderr)
            continue
        target[chr_name] = ids
    return target


def count_unitigs(unitigs: Path) -> int:
    n = 0
    with open_text(unitigs) as fin:
        for line in fin:
            if line.startswith(">"):
                n += 1
    return n


def parse_p_value(header_text: str) -> str:
    for part in header_text.split():
        if part.startswith("p="):
            return part[2:]
    return "NA"


def parse_male_female(header_text: str):
    male = None
    female = None
    for part in header_text.split():
        if part.startswith("male="):
            try:
                male = int(part.split("=", 1)[1])
            except ValueError:
                pass
        elif part.startswith("female="):
            try:
                female = int(part.split("=", 1)[1])
            except ValueError:
                pass
    return male, female


def longest_unitig_for_chr(intervals_by_chr, chr_name: str, lengths_map, headers_map, require_male_biased: bool = False):
    if chr_name not in intervals_by_chr:
        return None
    best_biased = None
    best_any = None
    for start, end, uid in intervals_by_chr[chr_name]:
        length = lengths_map.get(uid, end - start)
        header_text = headers_map.get(uid, uid)
        male, female = parse_male_female(header_text)
        candidate = (uid, length, start, end, header_text, male, female)
        if best_any is None or length > best_any[1]:
            best_any = candidate
        if male is not None and female is not None and male > female:
            if best_biased is None or length > best_biased[1]:
                best_biased = candidate
    chosen = best_biased if require_male_biased and best_biased else best_any
    if chosen is None:
        return None
    uid, length, start, end, header_text, male, female = chosen
    pval = parse_p_value(header_text)
    return uid, length, start, end, pval, header_text


def filter_male_biased(unitigs: Path, tmpdir: tempfile.TemporaryDirectory):
    out_path = Path(tmpdir.name) / f"{unitigs.name}.male_biased.fa"
    kept = 0
    skipped = 0
    with open_text(unitigs) as fin, out_path.open("w") as fout:
        header_text = None
        header_id = None
        seq_lines = []

        def flush():
            nonlocal kept, skipped
            if header_id is None:
                return
            seq = "".join(seq_lines).strip()
            if not seq:
                return
            male, female = parse_male_female(header_text)
            if male is not None and female is not None and male > female:
                fout.write(f">{header_text}\n{seq}\n")
                kept += 1
            else:
                skipped += 1

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header_text = line[1:].strip()
                header_id = header_text.split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        flush()
    print(f"[INFO] Filtered male-biased unitigs: kept {kept}, skipped {skipped}", file=sys.stderr)
    return out_path


def read_unitig_metadata(unitigs: Path):
    lengths = {}
    headers = {}
    total_bases = 0
    with open_text(unitigs) as fin:
        header_id = None
        header_text = None
        seq_lines = []

        def flush():
            nonlocal total_bases
            if header_id is None:
                return
            seq = "".join(seq_lines).strip()
            if not seq:
                return
            total_bases += len(seq)
            lengths[header_id] = len(seq)
            headers[header_id] = header_text

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header_text = line[1:].strip()
                header_id = header_text.split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        flush()
    return lengths, headers, total_bases


def write_unitigs_subset(unitigs: Path, ids: set, out_path: Path):
    if not ids:
        print(f"[WARN] No unitigs to write to {out_path}; skipping.", file=sys.stderr)
        return
    written = 0
    with open_text(unitigs) as fin, out_path.open("w") as fout:
        header_text = None
        header_id = None
        seq_lines = []
        def flush():
            nonlocal written
            if header_id is None or header_id not in ids:
                return
            seq = "".join(seq_lines).strip()
            if seq:
                fout.write(f">{header_text}\n{seq}\n")
                written += 1
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header_text = line[1:].strip()
                header_id = header_text.split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        flush()
    print(f"[INFO] Wrote {written} unitigs to {out_path}", file=sys.stderr)


def export_unitigs_to_fasta(unitigs: Path, chr_unitigs: OrderedDict, prefix: Path):
    """
    Write unitig sequences that mapped to selected chromosomes into separate FASTA files.
    """
    if not chr_unitigs:
        return
    union_ids = set()
    for ids in chr_unitigs.values():
        union_ids.update(ids)
    if not union_ids:
        print("[WARN] No unitigs matched requested chromosomes; FASTA export skipped.", file=sys.stderr)
        return

    out_handles = {}
    out_paths = {}
    counts = {chr_name: 0 for chr_name in chr_unitigs}
    try:
        for chr_name in chr_unitigs:
            out_path = prefix.parent / f"{prefix.name}.{chr_name}.fa"
            out_handles[chr_name] = out_path.open("w")
            out_paths[chr_name] = out_path
        with open_text(unitigs) as fin:
            header = None
            seq_lines = []

            def flush_sequence():
                if header is None or header not in union_ids:
                    return
                seq = "".join(seq_lines).strip()
                if not seq:
                    return
                for chr_name, ids in chr_unitigs.items():
                    if header in ids:
                        out_handles[chr_name].write(f">{header}\n{seq}\n")
                        counts[chr_name] += 1

            for line in fin:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    flush_sequence()
                    header = line[1:].strip().split()[0]
                    seq_lines = []
                else:
                    seq_lines.append(line.strip())
            flush_sequence()
    finally:
        for fh in out_handles.values():
            fh.close()
    for chr_name, out_path in out_paths.items():
        print(f"[INFO] Exported {counts[chr_name]} unitigs mapped to {chr_name} to {out_path}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(
        description="Map unitigs with minimap2 and/or BWA and compute windowed coverage."
    )
    ap.add_argument("--ref-a", dest="ref_a",
                    help="Reference genome FASTA for panel A (.fa, uncompressed for BWA)")
    ap.add_argument("--ref", dest="ref_a",
                    help=argparse.SUPPRESS)
    ap.add_argument("--unitigs", required=True,
                    help="Unitig FASTA (.fa[.gz])")
    ap.add_argument("--outdir", required=True,
                    help="Output directory for PNG/PDF and TSV")
    ap.add_argument("--window", type=int, default=5000,
                    help="Window size in bp (default: 5000)")
    ap.add_argument("--minimap2", default="minimap2",
                    help="Path to minimap2 binary")
    ap.add_argument("--bwa", default="bwa",
                    help="Path to bwa binary")
    ap.add_argument("--bwa-seed-length", type=int, default=31,
                    help="Seed length for bwa aln (-l, default: 31)")
    ap.add_argument("--threads", type=int, default=4,
                    help="Number of threads for mappers (used for all panels)")
    ap.add_argument("--mode", choices=["minimap", "minimap-strict", "bwa", "hybrid"], default="minimap",
                    help="Mapping mode (used for all panels): minimap, minimap-strict, bwa, or hybrid")
    ap.add_argument("--top-windows", type=int, default=0,
                    help="Report top N windows by coverage with unitig IDs (default: 0 = disable)")
    ap.add_argument("--export-chr", nargs="+", default=[],
                    help="Chromosome names for which mapped unitigs should be written to FASTA (one file per chromosome)")
    ap.add_argument("--ref-b",
                    help="Reference genome FASTA for panel B (.fa, uncompressed for BWA)")
    ap.add_argument("--unitigs-b",
                    help="Unitig FASTA for panel B (.fa[.gz])")
    ap.add_argument("--window-b", type=int,
                    help="Window size in bp for panel B (default: same as panel A)")
    ap.add_argument("--minimap2-b",
                    help="Path to minimap2 binary for panel B (default: same as panel A)")
    ap.add_argument("--bwa-b",
                    help="Path to bwa binary for panel B (default: same as panel A)")
    ap.add_argument("--bwa-seed-length-b", type=int,
                    help="Seed length for bwa aln in panel B (-l, default: panel A)")
    ap.add_argument("--top-windows-b", type=int,
                    help="Report top N windows for panel B (default: panel A)")
    ap.add_argument("--export-chr-b", nargs="+",
                    help="Chromosome names for FASTA export in panel B (default: panel A list)")
    ap.add_argument("--x-scale", type=float, default=0.85,
                    help="Scale factor for figure width (x-axis); <1 compresses, >1 stretches (default: 0.85)")
    ap.add_argument("--zoom-chr-a",
                    help="Chromosome name for an additional zoomed plot for panel A")
    ap.add_argument("--zoom-chr-b",
                    help="Chromosome name for an additional zoomed plot for panel B")
    ap.add_argument("--zoom-range-a",
                    help="Zoom range for panel A chromosome in Mb, e.g. 6-12")
    ap.add_argument("--zoom-range-b",
                    help="Zoom range for panel B chromosome in Mb, e.g. 2-8")
    ap.add_argument("--zoom-cut-a",
                    help="Cut range (Mb) to remove from panel A zoom/inset x-axis, e.g. 5-10")
    ap.add_argument("--zoom-cut-b",
                    help="Cut range (Mb) to remove from panel B zoom/inset x-axis, e.g. 200-400")
    ap.add_argument("--zoom-highlight-fill", action="store_true",
                    help="Highlight zoomed region on the main panel with a soft fill instead of a dashed box")
    ap.add_argument("--male-biased-only", action="store_true",
                    help="Filter input unitigs to male-biased entries (male > female in header) before mapping")
    ap.add_argument("--zoom-only", action="store_true",
                    help="Only render zoomed chromosome plots and skip the genome-wide plot")
    ap.add_argument("--insert-mode", action="store_true",
                    help="Embed zooms as inset plots in the main panels instead of separate columns")
    args = ap.parse_args()

    if not args.ref_a:
        print("[ERROR] Please provide --ref-a for panel A (alias --ref is deprecated).", file=sys.stderr)
        sys.exit(1)
    if args.x_scale <= 0:
        print("[ERROR] --x-scale must be > 0", file=sys.stderr)
        sys.exit(1)
    zoom_range_a = None
    zoom_range_b = None
    if args.zoom_range_a:
        try:
            zoom_range_a = parse_zoom_range(args.zoom_range_a)
        except ValueError as e:
            print(f"[ERROR] Invalid --zoom-range-a: {e}", file=sys.stderr)
            sys.exit(1)
    if args.zoom_range_b:
        try:
            zoom_range_b = parse_zoom_range(args.zoom_range_b)
        except ValueError as e:
            print(f"[ERROR] Invalid --zoom-range-b: {e}", file=sys.stderr)
            sys.exit(1)
    cut_range_a = None
    cut_range_b = None
    if args.zoom_cut_a:
        try:
            cut_range_a = parse_zoom_range(args.zoom_cut_a)
        except ValueError as e:
            print(f"[ERROR] Invalid --zoom-cut-a: {e}", file=sys.stderr)
            sys.exit(1)
    if args.zoom_cut_b:
        try:
            cut_range_b = parse_zoom_range(args.zoom_cut_b)
        except ValueError as e:
            print(f"[ERROR] Invalid --zoom-cut-b: {e}", file=sys.stderr)
            sys.exit(1)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    panels = [
        {
            "label": "A",
            "ref": Path(args.ref_a),
            "unitigs": Path(args.unitigs),
            "window": args.window,
            "minimap2": args.minimap2,
            "bwa": args.bwa,
            "threads": args.threads,
            "mode": args.mode,
            "top_windows": args.top_windows,
            "export_chr": args.export_chr,
            "bwa_seed_length": args.bwa_seed_length,
            "male_biased_only": args.male_biased_only,
        }
    ]

    if args.ref_b or args.unitigs_b:
        if not args.ref_b or not args.unitigs_b:
            print("[ERROR] Panel B requires both --ref-b and --unitigs-b", file=sys.stderr)
            sys.exit(1)
        panels.append(
            {
                "label": "B",
                "ref": Path(args.ref_b),
                "unitigs": Path(args.unitigs_b),
                "window": args.window_b if args.window_b is not None else args.window,
                "minimap2": args.minimap2_b if args.minimap2_b is not None else args.minimap2,
                "bwa": args.bwa_b if args.bwa_b is not None else args.bwa,
                "threads": args.threads,
                "mode": args.mode,
                "top_windows": args.top_windows_b if args.top_windows_b is not None else args.top_windows,
                "export_chr": args.export_chr_b if args.export_chr_b is not None else args.export_chr,
                "bwa_seed_length": args.bwa_seed_length_b if args.bwa_seed_length_b is not None else args.bwa_seed_length,
                "male_biased_only": args.male_biased_only,
            }
        )

    use_label_prefix = len(panels) > 1
    panel_results = []
    for panel in panels:
        panel_results.append(
            process_panel(
                panel["label"],
                panel["ref"],
                panel["unitigs"],
                outdir,
                panel["window"],
                panel["minimap2"],
                panel["bwa"],
                panel["threads"],
                panel["mode"],
                panel["top_windows"],
                panel["export_chr"],
                panel["bwa_seed_length"],
                panel["male_biased_only"],
                use_label_prefix,
            )
        )

    zoom_targets = {"A": "19", "B": "2"}
    if args.zoom_chr_a is not None:
        zoom_targets["A"] = args.zoom_chr_a
    if args.zoom_chr_b is not None:
        zoom_targets["B"] = args.zoom_chr_b
    zoom_ranges = {"A": zoom_range_a, "B": zoom_range_b}
    cut_ranges = {"A": cut_range_a, "B": cut_range_b}
    if len(panel_results) == 1:
        zoom_targets["B"] = None
        zoom_ranges["B"] = None
        cut_ranges["B"] = None

    if not args.zoom_only:
        if any(zoom_targets.values()):
            if len(panel_results) == 1:
                plot_prefix = panel_results[0]["prefix"]
            else:
                plot_prefix = outdir / f"{panel_results[0]['prefix'].name}__{panel_results[1]['prefix'].name}"
            out_png = plot_prefix.with_suffix(".png")
            out_pdf = plot_prefix.with_suffix(".pdf")
            print(f"[INFO] Plot PNG (with zooms): {out_png}", file=sys.stderr)
            print(f"[INFO] Plot PDF (with zooms): {out_pdf}", file=sys.stderr)
            plot_coverage_with_zooms(
                panel_results, zoom_targets, zoom_ranges, cut_ranges, out_png, out_pdf, x_scale=args.x_scale, insert_mode=args.insert_mode, zoom_fill=args.zoom_highlight_fill
            )
        else:
            if len(panel_results) == 1:
                plot_prefix = panel_results[0]["prefix"]
            else:
                plot_prefix = outdir / f"{panel_results[0]['prefix'].name}__{panel_results[1]['prefix'].name}"
            out_png = plot_prefix.with_suffix(".png")
            out_pdf = plot_prefix.with_suffix(".pdf")
            print(f"[INFO] Plot PNG: {out_png}", file=sys.stderr)
            print(f"[INFO] Plot PDF: {out_pdf}", file=sys.stderr)
            plot_coverage_panels(panel_results, out_png, out_pdf, x_scale=args.x_scale)

    for panel in panel_results:
        chrom_target = zoom_targets.get(panel["label"])
        if not chrom_target:
            continue
        out_zoom_png, out_zoom_pdf = make_zoom_paths(panel["prefix"], chrom_target)
        plot_zoom_chromosome(
            panel,
            chrom_target,
            out_zoom_png,
            out_zoom_pdf,
            zoom_range_mb=zoom_ranges.get(panel["label"]),
            cut_range_mb=cut_ranges.get(panel["label"]) if cut_ranges else None,
        )


if __name__ == "__main__":
    main()
