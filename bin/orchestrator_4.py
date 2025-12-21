#!/usr/bin/env python3
import argparse
import math
import gzip
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path
from datetime import datetime
from typing import Dict, Iterable, List

def ts():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def make_logger(log_unitig_path: Path):
    log_unitig_path.parent.mkdir(parents=True, exist_ok=True)

    def log(msg: str, to_main: bool = False):
        line = f"[{ts()}] [INFO] {msg}"
        print(line, flush=True)
        with log_unitig_path.open("a", encoding="utf-8") as fh:
            fh.write(line + "\n")
    return log

def read_present_count(bin_path: Path) -> int:
    if not bin_path.is_file():
        return 0
    try:
        with bin_path.open("rb") as fh:
            header = fh.read(4)
            if len(header) != 4:
                return 0
            return int.from_bytes(header, "little")
    except OSError:
        return 0

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--metadata", default="metadata.tsv")
    p.add_argument("--outdir", default=None)
    p.add_argument("--threads", type=int, default=3)
    p.add_argument("--carrier-threshold", type=int, default=3)
    p.add_argument("--min-carrier", dest="carrier_threshold", type=int, help=argparse.SUPPRESS)
    p.add_argument("--anchor-bp", type=int, default=12)
    p.add_argument("--query-split-size", type=int, default=32767)
    p.add_argument("--kmer-length", type=int, default=19)
    p.add_argument("--min-mem-length", type=int, default=None,
                   help="Minimum MEM length (default: k+1 from unitig stats or kmer-length+1)")
    p.add_argument("--reverse-complement", default="for_reference")
    p.add_argument("--no-index-cache", action="store_true",
                   help="Disable nt_mini_matcher index cache (default: build once and reuse)")
    p.add_argument("--augment-k-singletons", action="store_true",
                   help="Append k-length unitigs from candidate_kmers to presence table (default: off)")
    p.add_argument("--short-frac", type=float, default=1.0)
    p.add_argument("--medium-frac", type=float, default=0.95)
    p.add_argument("--long-frac", type=float, default=0.9)
    p.add_argument("--short-min-bp", type=int, default=None)
    p.add_argument("--medium-min-bp", type=int, default=None)
    p.add_argument("--long-min-bp", type=int, default=None)
    p.add_argument("--bridge-min-span", type=int, default=150)
    p.add_argument("--short-min-segment", type=int, default=None)
    p.add_argument("--medium-min-segment", type=int, default=None)
    p.add_argument("--long-min-segment", type=int, default=None)
    p.add_argument("--min-bridging-pairs", type=int, default=2)
    p.add_argument("--full-gapless-long", action="store_true")
    p.add_argument("--max-replicate-ref", type=int, default=None,
                   help="Limit replicated seeds on reference (nt_mini_matcher --max_replicate_ref)")
    p.add_argument("--max-replicate-qry", type=int, default=None,
                   help="Limit replicated seeds on query (nt_mini_matcher --max_replicate_qry)")
    return p.parse_args()

def load_metadata(meta_path: Path):
    samples = {}
    traits = set()
    with meta_path.open("r", encoding="utf-8") as fh:
        header = fh.readline()
        cols = [c.strip() for c in header.rstrip("\n").split("\t")]
        idx_path = cols.index("path")
        idx_sample = cols.index("sample_id")
        idx_trait = cols.index("trait")
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx_path, idx_sample, idx_trait):
                continue
            p = parts[idx_path].strip()
            sid = parts[idx_sample].strip()
            tr = parts[idx_trait].strip()
            if not sid or not p or not tr:
                continue
            samples.setdefault(sid, {"trait": tr, "paths": []})
            samples[sid]["paths"].append(p)
            traits.add(tr)
    return samples, sorted(traits)

def parse_unitig_stats(stats_path: Path):
    if not stats_path.is_file():
        return None
    total = None
    median = None
    max_len = None
    k_len = None
    minlen_threshold = None
    singleton_k = None
    with stats_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            line = line.strip()
            if not line or "\t" not in line:
                continue
            key, value = line.split("\t", 1)
            if key.startswith("n_unitigs_after_minlen_"):
                try:
                    total = int(value)
                except ValueError:
                    total = None
                try:
                    minlen_threshold = int(key.rsplit("_", 1)[1])
                except (IndexError, ValueError):
                    minlen_threshold = None
            elif key in {"n_unitigs_for_matching", "n_unitigs"}:
                try:
                    total = int(value)
                except ValueError:
                    total = None
            elif key == "median_len":
                try:
                    median = float(value)
                except ValueError:
                    median = None
            elif key == "max_len":
                try:
                    max_len = int(value)
                except ValueError:
                    max_len = None
            elif key == "k":
                try:
                    k_len = int(value)
                except ValueError:
                    k_len = None
            elif key in {"n_singleton_k", "unitigs_length_k", "unitigs_of_length_k"}:
                try:
                    singleton_k = int(value)
                except ValueError:
                    singleton_k = None
    return {
        "n_unitigs": total,
        "median_len": median,
        "max_len": max_len,
        "k_len": k_len,
        "n_singleton_k": singleton_k,
        "minlen_threshold": minlen_threshold,
    }

def load_length_histogram(stats_path: Path):
    hist = {}
    if not stats_path.is_file():
        return hist
    in_hist = False
    with stats_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("# length_histogram"):
                in_hist = True
                continue
            if not in_hist:
                continue
            if line.startswith("#") or line.startswith("length"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                continue
            try:
                L = int(parts[0]); c = int(parts[1])
            except ValueError:
                continue
            if L < 0 or c < 0:
                continue
            hist[L] = hist.get(L, 0) + c
    return hist

def weighted_quantile_from_hist(hist: dict, q: float):
    if not hist:
        return None
    total = sum(hist.values())
    if total <= 0:
        return None
    target = q * total
    acc = 0
    for L in sorted(hist.keys()):
        acc += hist[L]
        if acc >= target:
            return L
    return max(hist.keys())

def weighted_median_from_hist(hist: dict):
    return weighted_quantile_from_hist(hist, 0.5)

def compute_class_median(hist: dict, lo: int, hi: int):
    sub = {L: c for L, c in hist.items() if lo <= L <= hi}
    return weighted_median_from_hist(sub)

def compute_length_params_heuristic(hist: dict, k_len: int, fallback_median: float | None, fallback_max: int | None):
    legacy_short_max = int(round(fallback_median)) if fallback_median else 40
    legacy_medium_max = int(round((fallback_median or 40) * 2.5))
    if fallback_max and legacy_medium_max > fallback_max:
        legacy_medium_max = fallback_max
    legacy_short_min_bp, legacy_med_min_bp, legacy_long_min_bp = 37, 92, 150
    legacy_short_min_seg, legacy_med_min_seg, legacy_long_min_seg = 31, 50, 60

    if not hist or sum(hist.values()) <= 0:
        return (
            legacy_short_max,
            legacy_medium_max,
            legacy_short_min_bp,
            legacy_med_min_bp,
            legacy_long_min_bp,
            legacy_short_min_seg,
            legacy_med_min_seg,
            legacy_long_min_seg,
        )

    q25 = weighted_quantile_from_hist(hist, 0.25) or legacy_short_max
    q75 = weighted_quantile_from_hist(hist, 0.75) or legacy_medium_max

    short_max = max(k_len + 1, q25)
    medium_max = max(short_max + 1, q75)

    L_short_med = compute_class_median(hist, min(hist.keys()), short_max) or short_max
    L_med_med = compute_class_median(hist, short_max + 1, medium_max) or medium_max
    L_long_med = compute_class_median(hist, medium_max + 1, max(hist.keys())) or medium_max

    short_min_bp = max(k_len, math.ceil(0.9 * L_short_med))
    med_min_bp = max(short_min_bp + 1, math.ceil(0.8 * L_med_med))
    long_min_bp = max(med_min_bp + 1, math.ceil(0.8 * L_long_med))

    short_min_seg = short_min_bp
    med_min_seg = max(short_min_bp, min(med_min_bp, L_med_med))
    long_min_seg = max(med_min_bp, min(long_min_bp, L_long_med))

    return (
        short_max,
        medium_max,
        short_min_bp,
        med_min_bp,
        long_min_bp,
        short_min_seg,
        med_min_seg,
        long_min_seg,
    )


def sanitize_trait(s: str) -> str:
    return "".join(ch for ch in s if ch.isalnum())

def unique_prefixes(trait0: str, trait1: str):
    a = sanitize_trait(trait0)
    b = sanitize_trait(trait1)
    if not a and not b:
        return "0", "1"
    if not a:
        a = "0"
    if not b:
        b = "1"
    al = a.lower()
    bl = b.lower()
    L = max(len(al), len(bl))
    for i in range(1, L + 1):
        pa = al[:i]
        pb = bl[:i]
        if pa != pb:
            return a[: min(i, len(a))], b[: min(i, len(b))]
    return a + "0", b + "1"

def ordered_traits_from_counts(trait_counts: Counter):
    ordered = sorted(trait_counts.items(), key=lambda kv: (-kv[1], kv[0]))
    if len(ordered) < 2:
        raise SystemExit("Expected exactly two traits in metadata.")
    return ordered[0][0], ordered[1][0]

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

def lengths_from_presence(presence_path: Path, len_map: Dict[str, int]) -> List[int]:
    lengths: List[int] = []
    if not presence_path.is_file():
        return lengths
    with presence_path.open("r", encoding="utf-8") as fh:
        _ = fh.readline()
        for line in fh:
            if not line.strip():
                continue
            uid = line.split("\t", 1)[0]
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

def parse_candidate_header(h: str, tag0: str, tag1: str):
    h = h.lstrip(">")
    m = re.match(r"([0-9]+)(.*)", h)
    if not m:
        return None
    rest = m.group(2)
    if not rest.startswith(tag0):
        return None
    rest = rest[len(tag0) :]
    m = re.match(r"([0-9]+)(.*)", rest)
    if not m:
        return None
    prev0 = int(m.group(1))
    rest = m.group(2)
    if not rest.startswith(tag1):
        return None
    prev1_part = rest[len(tag1) :]
    if not prev1_part.isdigit():
        return None
    prev1 = int(prev1_part)
    return prev0, prev1

def fasta_records(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    header = None
    seq_lines = []
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header.strip(), "".join(seq_lines).replace(" ", "").replace("\t", "").strip()
                header = line.strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            yield header.strip(), "".join(seq_lines).replace(" ", "").replace("\t", "").strip()

def build_singleton_presence(unitigs_fa: Path,
                             candidate_fa: Path,
                             k_len: int,
                             trait0: str,
                             trait1: str,
                             traits_header,
                             log):
    if k_len is None:
        log("[SINGLETON] Skipping augmentation (k not found in unitig stats)", to_main=True)
        return []

    tag0, tag1 = unique_prefixes(trait0, trait1)

    seq_to_prev = {}
    for h, seq in fasta_records(candidate_fa):
        if len(seq) != k_len:
            continue
        parsed = parse_candidate_header(h, tag0, tag1)
        if parsed is None:
            continue
        prev0, prev1 = parsed
        seq_to_prev[seq] = {trait0: prev0, trait1: prev1}

    if not seq_to_prev:
        log("[SINGLETON] No k-length candidates found to augment.", to_main=True)
        return []

    missing_rows = []
    for h, seq in fasta_records(unitigs_fa):
        if len(seq) != k_len:
            continue
        unitig_id = h.lstrip(">").split()[0]
        prev = seq_to_prev.get(seq)
        if prev is None:
            continue
        row = [unitig_id]
        for t in traits_header:
            row.append(str(prev.get(t, 0)))
        missing_rows.append("\t".join(row))

    return missing_rows

def read_summary_stats(stats_path: Path):
    if not stats_path.is_file():
        return {}
    res = {}
    trait_bias = {}
    full_prev = {}
    trait_specific_only = {}
    full_prev_absent_others = {}
    with stats_path.open("r", encoding="utf-8") as fh:
        mode_presence_hist = False
        mode_absdiff_hist = False
        presence_bins = set()
        absdiff_bins = set()
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            if line.startswith("total_unitigs"):
                _, v = line.split("\t")
                res["total_unitigs"] = int(v)
            elif line.startswith("present_unitigs"):
                _, v = line.split("\t")
                res["present_unitigs"] = int(v)
            elif line.startswith("absent_unitigs"):
                _, v = line.split("\t")
                res["absent_unitigs"] = int(v)
            elif line == "presence_total_hist":
                mode_presence_hist = True
                mode_absdiff_hist = False
            elif line == "abs_diff_hist":
                mode_presence_hist = False
                mode_absdiff_hist = True
            elif mode_presence_hist:
                if line.startswith("total_presence"):
                    continue
                parts = line.split("\t")
                if len(parts) == 2:
                    tp = int(parts[0])
                    presence_bins.add(tp)
            elif mode_absdiff_hist:
                if line.startswith("abs_diff"):
                    continue
                parts = line.split("\t")
                if len(parts) == 2:
                    d = int(parts[0])
                    absdiff_bins.add(d)
            elif line.startswith("TRAIT_BIAS"):
                parts = line.split("\t")
                if len(parts) >= 3:
                    trait = parts[1]
                    rest = parts[2]
                    vals = {}
                    for token in rest.split():
                        if "=" in token:
                            k, v = token.split("=", 1)
                            vals[k] = int(v)
                    trait_bias[trait] = vals
            elif line.startswith("FULL_PREVALENCE_ABSENT_OTHERS"):
                _, trait, v = line.split("\t")
                full_prev_absent_others[trait] = int(v)
            elif line.startswith("FULL_PREVALENCE\t"):
                _, trait, v = line.split("\t")
                full_prev[trait] = int(v)
            elif line.startswith("TRAIT_SPECIFIC_ONLY"):
                _, trait, v = line.split("\t")
                trait_specific_only[trait] = int(v)
    if presence_bins:
        res["presence_total_hist_bins"] = len(presence_bins)
        res["presence_total_hist_max_total_presence"] = max(presence_bins)
    if absdiff_bins:
        res["abs_diff_hist_bins"] = len(absdiff_bins)
        res["abs_diff_hist_max_abs_diff"] = max(absdiff_bins)
    if trait_bias:
        res["TRAIT_BIAS"] = trait_bias
    if full_prev:
        res["FULL_PREVALENCE"] = full_prev
    if trait_specific_only:
        res["TRAIT_SPECIFIC_ONLY"] = trait_specific_only
    if full_prev_absent_others:
        res["FULL_PREVALENCE_ABSENT_OTHERS"] = full_prev_absent_others
    return res

def main():
    args = parse_args()
    root_dir = Path(__file__).resolve().parent.parent

    if args.outdir is None:
        outdir = root_dir
    else:
        outdir = root_dir / args.outdir if not Path(args.outdir).is_absolute() else Path(args.outdir)

    results_dir = outdir / "results"
    logging_dir = outdir / "logging"
    stats_dir = results_dir / "statistics"
    matcher_dir = results_dir / "unitig_matcher"

    results_dir.mkdir(parents=True, exist_ok=True)
    logging_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)
    matcher_dir.mkdir(parents=True, exist_ok=True)

    log_unitig_path = logging_dir / "unitig_matching.log"
    log = make_logger(log_unitig_path)

    meta_path = root_dir / args.metadata
    if not meta_path.is_file():
        raise SystemExit(f"Metadata not found: {meta_path}")

    unitigs_fa = results_dir / "unitigs.fa"
    if not unitigs_fa.is_file():
        raise SystemExit(f"unitigs.fa not found at {unitigs_fa}")

    unitig_stats_path = stats_dir / "unitig_statistics.txt"
    stats = parse_unitig_stats(unitig_stats_path) if unitig_stats_path.is_file() else {}
    length_hist = load_length_histogram(unitig_stats_path) if unitig_stats_path.is_file() else {}
    k_len = stats.get("k_len") if stats else None
    (
        short_max_len,
        medium_max_len,
        heur_short_min_bp,
        heur_med_min_bp,
        heur_long_min_bp,
        heur_short_min_seg,
        heur_med_min_seg,
        heur_long_min_seg,
    ) = compute_length_params_heuristic(
        length_hist,
        k_len or args.kmer_length,
        stats.get("median_len") if stats else None,
        stats.get("max_len") if stats else None,
    )

    threads_eff = min(args.threads, 16)

    if args.short_min_bp is None:
        args.short_min_bp = heur_short_min_bp
    if args.medium_min_bp is None:
        args.medium_min_bp = heur_med_min_bp
    if args.long_min_bp is None:
        args.long_min_bp = heur_long_min_bp
    if args.short_min_segment is None:
        args.short_min_segment = heur_short_min_seg
    if args.medium_min_segment is None:
        args.medium_min_segment = heur_med_min_seg
    if args.long_min_segment is None:
        args.long_min_segment = heur_long_min_seg

    if args.min_mem_length is None:
        args.min_mem_length = (k_len or args.kmer_length) + 1
        log(
            f"min_mem_length not set; using k+1 = {args.min_mem_length} "
            f"(k from stats={k_len or 'n/a'}, kmer_length arg={args.kmer_length})",
            to_main=True,
        )

    short_frac = args.short_frac
    medium_frac = args.medium_frac
    long_frac = args.long_frac

    samples, traits = load_metadata(meta_path)
    n_samples = len(samples)
    trait_counts = {}
    all_even = True
    for sid, info in samples.items():
        n = len(info["paths"])
        if n % 2 != 0:
            all_even = False
        trait_counts[info["trait"]] = trait_counts.get(info["trait"], 0) + 1
    trait0, trait1 = ordered_traits_from_counts(Counter(trait_counts))

    if all_even:
        mode_global = "paired-end (auto, all samples have even #FASTQ)"
    else:
        mode_global = "mixed"

    log("=== UNITIG MATCHING START ===", to_main=True)
    log(f"METADATA\t{meta_path.name}", to_main=True)
    log("=== PARAMETERS ===", to_main=True)
    log(f"kmer_length\t{args.kmer_length}", to_main=True)
    log(f"min_mem_length\t{args.min_mem_length}", to_main=True)
    log(f"threads\t{threads_eff}", to_main=True)
    log(f"reverse_complement\t{args.reverse_complement}", to_main=True)
    log(f"carrier_threshold\t{args.carrier_threshold}", to_main=True)
    log(f"augment_k_singletons\t{args.augment_k_singletons}", to_main=True)
    log(f"short_max\t{short_max_len}", to_main=True)
    log(f"medium_max\t{medium_max_len}", to_main=True)
    log(f"short_frac\t{short_frac}", to_main=True)
    log(f"medium_frac\t{medium_frac}", to_main=True)
    log(f"long_frac\t{long_frac}", to_main=True)
    log(f"short_min_bp\t{args.short_min_bp}", to_main=True)
    log(f"medium_min_bp\t{args.medium_min_bp}", to_main=True)
    log(f"long_min_bp\t{args.long_min_bp}", to_main=True)
    log(f"bridge_min_span\t{args.bridge_min_span}", to_main=True)
    log(f"short_min_segment\t{args.short_min_segment}", to_main=True)
    log(f"medium_min_segment\t{args.medium_min_segment}", to_main=True)
    log(f"long_min_segment\t{args.long_min_segment}", to_main=True)
    log(f"min_bridging_pairs\t{args.min_bridging_pairs}", to_main=True)
    log(f"full_gapless_long\t{args.full_gapless_long}", to_main=True)
    log(f"query_split_size\t{args.query_split_size}", to_main=True)
    log(f"index_cache\t{'disabled' if args.no_index_cache else 'auto'}", to_main=True)
    if args.max_replicate_ref is not None:
        log(f"max_replicate_ref\t{args.max_replicate_ref}", to_main=True)
    if args.max_replicate_qry is not None:
        log(f"max_replicate_qry\t{args.max_replicate_qry}", to_main=True)
    # mode_global logging removed per request

    log("=== UNITIG LENGTH CLASSES ===", to_main=True)
    log(f"short\tlen <= {short_max_len}", to_main=True)
    log(f"medium\t{short_max_len+1} <= len <= {medium_max_len}", to_main=True)
    log(f"long\tlen > {medium_max_len}", to_main=True)

    log(f"n_samples\t{n_samples}", to_main=True)
    for t in traits:
        log(f"trait_samples\t{t}\t{trait_counts.get(t, 0)}", to_main=True)

    index_cache_path = matcher_dir / "unitigs.idx"

    use_cache = not args.no_index_cache
    if args.no_index_cache:
        log("[CACHE] disabled (no-index-cache)", to_main=True)
    else:
        build_cmd = [
            str(root_dir / "bin" / "nt_mini_matcher.x"),
            "--index-cache", str(index_cache_path),
            "--build-index-only",
            "--kmer_length", str(args.kmer_length),
            "--minimum_mem_length", str(args.min_mem_length),
            "--reverse_complement", args.reverse_complement,
            str(unitigs_fa),
        ]
        log(f"[CACHE] build index -> {index_cache_path.name}")
        res_build = subprocess.run(build_cmd, cwd=matcher_dir,
                                   capture_output=True, text=True)
        if res_build.returncode != 0:
            err_out = res_build.stdout or ""
            err_err = res_build.stderr or ""
            log(
                f"[CACHE] build failed rc={res_build.returncode} stdout={err_out!r} stderr={err_err!r}",
                to_main=True,
            )
            raise SystemExit(
                f"nt_mini_matcher cache build failed with exit code {res_build.returncode}"
            )
        if not index_cache_path.is_file():
            raise SystemExit(
                f"nt_mini_matcher cache build reported success but file missing: {index_cache_path}"
            )

    log("=== PROCESSING SAMPLES ===", to_main=True)

    orchestrator_root = root_dir
    sorted_ids = sorted(samples.keys())

    for sid in sorted_ids:
        info = samples[sid]
        trait = info["trait"]
        fastqs = info["paths"]
        n_fastq = len(fastqs)

        log(f"[{sid}] START", to_main=True)

        if n_fastq == 0:
            log(f"[{sid}] SKIP (no FASTQs)", to_main=True)
            continue

        present_bin = matcher_dir / f"{sid}.present_unitigs.bin"

        paired = (n_fastq % 2 == 0 and n_fastq >= 2)

        base_cmd_post = [
            str(root_dir / "bin" / "postprocess_minimatcher"),
            str(unitigs_fa),
            "-",
            str(present_bin),
            "--carrier-threshold", str(args.carrier_threshold),
            "--short-max", str(short_max_len),
            "--medium-max", str(medium_max_len),
            "--short-frac", str(short_frac),
            "--medium-frac", str(medium_frac),
            "--long-frac", str(long_frac),
            "--anchor-bp", str(args.anchor_bp),
            "--short-min-bp", str(args.short_min_bp),
            "--medium-min-bp", str(args.medium_min_bp),
            "--long-min-bp", str(args.long_min_bp),
            "--bridge-min-span", str(args.bridge_min_span),
            "--short-min-segment", str(args.short_min_segment),
            "--medium-min-segment", str(args.medium_min_segment),
            "--long-min-segment", str(args.long_min_segment),
            "--min-bridging-pairs", str(args.min_bridging_pairs),
        ]
        if args.full_gapless_long:
            base_cmd_post.append("--full-gapless-long")

        base_cmd_nt = [
            str(root_dir / "bin" / "nt_mini_matcher.x"),
            "--threads", str(threads_eff),
            "--query_split_size", str(args.query_split_size),
            "--reverse_complement", args.reverse_complement,
            "--display", "q.seqlen q.read_pairs",
            "--kmer_length", str(args.kmer_length),
            "--minimum_mem_length", str(args.min_mem_length),
        ]
        if args.max_replicate_ref is not None:
            base_cmd_nt += ["--max_replicate_ref", str(args.max_replicate_ref)]
        if args.max_replicate_qry is not None:
            base_cmd_nt += ["--max_replicate_qry", str(args.max_replicate_qry)]
        if use_cache:
            base_cmd_nt += ["--index-cache", str(index_cache_path)]
        if paired:
            base_cmd_nt.append("--read_pairs")
        base_cmd_nt += [str(unitigs_fa)] + fastqs

        log(f"[{sid}] RUN nt_mini_matcher n_fastq={n_fastq} trait={trait}")
        log(f"[{sid}] RUN postprocessing trait={trait}")
        p_nt = subprocess.Popen(
            base_cmd_nt,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=matcher_dir,
            bufsize=0,
        )
        p_post = subprocess.Popen(base_cmd_post, stdin=p_nt.stdout, cwd=matcher_dir)
        if p_nt.stdout is not None:
            p_nt.stdout.close()

        post_rc = p_post.wait()
        stderr_nt = ""
        if p_nt.stderr is not None:
            try:
                stderr_nt = p_nt.stderr.read().decode(errors="ignore").strip()
            except Exception:
                stderr_nt = ""
        nt_rc = p_nt.wait()

        if nt_rc != 0:
            msg = f"nt_mini_matcher failed for {sid} with exit code {nt_rc}"
            if stderr_nt:
                msg += f" (stderr: {stderr_nt})"
            raise SystemExit(msg)
        if post_rc != 0:
            raise SystemExit(
                f"postprocess_minimatcher failed for {sid} "
                f"with exit code {post_rc}"
            )
        count_present = read_present_count(present_bin)
        log(
            f"[{sid}] PRESENT_UNITIGS\t{count_present}",
            to_main=True,
        )
        log(f"[{sid}] DONE", to_main=True)

    presence_table = results_dir / "unitig_presence.tsv"
    matches_stats = stats_dir / "unitig_matches.txt"

    log("[MERGE] START", to_main=True)

    cmd_merge = [
        str(orchestrator_root / "bin" / "merge_matches"),
        "--indir", str(matcher_dir),
        "--metadata", str(meta_path),
        "--out", str(presence_table),
        "--stats", str(matches_stats),
        "--unitig-stats", str(unitig_stats_path),
    ]

    subprocess.run(cmd_merge, check=True)

    log("[MERGE] DONE", to_main=True)

    if args.augment_k_singletons:
        if k_len is None:
            log("[SINGLETON] Skipping augmentation (k not found in unitig stats)", to_main=True)
        elif args.min_mem_length <= k_len:
            log(f"[SINGLETON] Skipping augmentation (min_mem_length <= k={k_len}, already matched via nt_mini_matcher)", to_main=True)
        elif args.min_mem_length > k_len + 1:
            log(f"[SINGLETON] Skipping augmentation (min_mem_length={args.min_mem_length} > k+1={k_len + 1})", to_main=True)
        else:
            candidate_fa = results_dir / "candidate_kmers.fasta"
            if not candidate_fa.is_file():
                cand_gz = candidate_fa.with_suffix(candidate_fa.suffix + ".gz")
                candidate_fa = cand_gz if cand_gz.is_file() else None

            if candidate_fa is None:
                log("[SINGLETON] candidate_kmers FASTA not found; skipping augmentation", to_main=True)
            else:
                traits_header = []
                existing_ids = set()
                if presence_table.is_file():
                    with presence_table.open("r", encoding="utf-8") as fh:
                        header_line = fh.readline().strip()
                        traits_header = header_line.split("\t")[1:] if header_line else []
                        for line in fh:
                            line = line.strip()
                            if not line:
                                continue
                            existing_ids.add(line.split("\t", 1)[0])
                if not traits_header:
                    traits_header = sorted(trait_counts.keys())
                    with presence_table.open("w", encoding="utf-8") as fh:
                        fh.write("unitig\t" + "\t".join(traits_header) + "\n")

                singleton_rows = build_singleton_presence(
                    unitigs_fa,
                    candidate_fa,
                    k_len,
                    trait0,
                    trait1,
                    traits_header,
                    log,
                )
                singleton_rows = [
                    r for r in singleton_rows
                    if r.split("\t", 1)[0] not in existing_ids
                ]
                if singleton_rows:
                    with presence_table.open("a", encoding="utf-8") as fh:
                        for row in singleton_rows:
                            fh.write(row + "\n")
                    log(f"[SINGLETON] Added {len(singleton_rows)} k-length unitigs from candidate_kmers", to_main=True)
                else:
                    log("[SINGLETON] No additional k-length unitigs added (already present or missing matches)", to_main=True)
    else:
        log("[SINGLETON] Augmentation disabled (default)", to_main=True)

    try:
        len_map = load_unitig_lengths(unitigs_fa)
        matched_lengths = lengths_from_presence(presence_table, len_map)
        append_histogram(matched_lengths, matches_stats, "present unitigs")
    except Exception as exc:
        log(f"[WARN] Failed to append present unitig length histogram: {exc}", to_main=True)

    try:
        summary = read_summary_stats(matches_stats)
        log("=== FINAL SUMMARY ===", to_main=True)
        if "total_unitigs" in summary:
            log(f"total_unitigs\t{summary['total_unitigs']}", to_main=True)
        if "present_unitigs" in summary:
            log(f"present_unitigs\t{summary['present_unitigs']}", to_main=True)
        if "absent_unitigs" in summary:
            log(f"absent_unitigs\t{summary['absent_unitigs']}", to_main=True)
        if "presence_total_hist_bins" in summary:
            log(f"presence_total_hist_bins\t{summary['presence_total_hist_bins']}", to_main=True)
        if "presence_total_hist_max_total_presence" in summary:
            log(f"presence_total_hist_max_total_presence\t{summary['presence_total_hist_max_total_presence']}", to_main=True)
        if "abs_diff_hist_bins" in summary:
            log(f"abs_diff_hist_bins\t{summary['abs_diff_hist_bins']}", to_main=True)
        if "abs_diff_hist_max_abs_diff" in summary:
            log(f"abs_diff_hist_max_abs_diff\t{summary['abs_diff_hist_max_abs_diff']}", to_main=True)

        trait_bias = summary.get("TRAIT_BIAS", {})
        for t, vals in trait_bias.items():
            specific_only = vals.get("specific_only", 0)
            biased_gt = vals.get("biased_gt_other", 0)
            log(f"TRAIT_BIAS\t{t}\tspecific_only={specific_only}\tbiased_gt_other={biased_gt}", to_main=True)

        full_prev = summary.get("FULL_PREVALENCE", {})
        for t, v in full_prev.items():
            log(f"FULL_PREVALENCE\t{t}\t{v}", to_main=True)

        trait_spec_only = summary.get("TRAIT_SPECIFIC_ONLY", {})
        for t, v in trait_spec_only.items():
            log(f"TRAIT_SPECIFIC_ONLY\t{t}\t{v}", to_main=True)

        full_prev_absent = summary.get("FULL_PREVALENCE_ABSENT_OTHERS", {})
        for t, v in full_prev_absent.items():
            log(f"FULL_PREVALENCE_ABSENT_OTHERS\t{t}\t{v}", to_main=True)
    except Exception as exc:
        log(f"[WARN] Failed to read/write summary stats: {exc}", to_main=True)

    log(f"[CLEANUP] Removing temporary matcher directory: {matcher_dir.name}", to_main=True)
    try:
        shutil.rmtree(matcher_dir)
        log(f"[CLEANUP] Removed {matcher_dir.name}", to_main=True)
    except FileNotFoundError:
        log(f"[CLEANUP] Skipped removal, not found: {matcher_dir.name}", to_main=True)
    except Exception as exc:
        log(f"[CLEANUP] Failed to remove {matcher_dir}: {exc}", to_main=True)

if __name__ == "__main__":
    main()
