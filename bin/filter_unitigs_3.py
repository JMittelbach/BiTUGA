#!/usr/bin/env python3
import argparse, gzip, logging, sys, time, os, shutil
from collections import Counter
from statistics import mean, median

LC_CHECK_HOMO = False
LC_CHECK_DI = False
LC_CHECK_TRI = False
from pathlib import Path

def ensure_dir(p: Path): p.mkdir(parents=True, exist_ok=True)

def open_text_read(p: Path):
    if str(p).endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p, "rt")

def open_text_write(p: Path):
    if str(p).endswith(".gz"):
        return gzip.open(p, "wt")
    return open(p, "wt")

def setup_logger(log_dir: Path):
    ensure_dir(log_dir)
    logger = logging.getLogger("filter_unitigs")
    if logger.handlers:
        logger.handlers.clear()
    logger.setLevel(logging.INFO)
    logger.propagate = False
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh = logging.FileHandler(log_dir / "build_unitigs.log", mode="a")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    return logger

def ts_line(msg: str) -> str:
    return f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] [INFO] {msg}"

def iter_fasta(path: Path):
    with open_text_read(path) as fh:
        h = None
        seq = []
        for line in fh:
            if not line:
                continue
            if line[0] == ">":
                if h is not None:
                    yield h, "".join(seq)
                h = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if h is not None:
            yield h, "".join(seq)

def n50_l50(lengths):
    if not lengths:
        return 0, 0
    lengths_sorted = sorted(lengths, reverse=True)
    half = sum(lengths_sorted) / 2.0
    acc = 0
    for idx, L in enumerate(lengths_sorted, start=1):
        acc += L
        if acc >= half:
            return L, idx
    return lengths_sorted[-1], len(lengths_sorted)

def lc_reason(seq: str):
    """Return reason string if sequence fails LC filters, else None."""
    s = seq.upper()
    L = len(s)
    if LC_CHECK_HOMO and L > 0 and len(set(s)) == 1:
        return "homopolymer"
    if LC_CHECK_DI and L >= 4:
        pat = s[:2]
        if pat * (L // 2) == s[: (L // 2) * 2] and L % 2 == 0:
            return "dinucleotide_repeat"
    if LC_CHECK_TRI and L >= 6:
        pat = s[:3]
        if pat * (L // 3) == s[: (L // 3) * 3] and L % 3 == 0:
            return "trinucleotide_repeat"
    return None

def main():
    ap = argparse.ArgumentParser(description="Filter unitigs by minimum length and low-complexity; overwrite unitigs.fa only if something was removed.")
    ap.add_argument("--base-dir", required=True)
    ap.add_argument("--min-length", type=int, required=True)
    args = ap.parse_args()

    base = Path(args.base_dir).resolve()
    res_dir = base / "results"
    log_dir = base / "logging"
    stats_dir = res_dir / "statistics"
    ensure_dir(res_dir); ensure_dir(log_dir); ensure_dir(stats_dir)

    logger = setup_logger(log_dir)
    main_log = log_dir / "main.log"

    in_fa = res_dir / "unitigs.fa"
    if not in_fa.exists():
        alt = res_dir / "unitigs.fasta"
        if alt.exists():
            in_fa = alt
        else:
            print(f"[filter_unitigs] input not found: {res_dir}/unitigs.fa", file=sys.stderr)
            sys.exit(2)

    thr = int(args.min_length)
    stats_out = stats_dir / "unitig_statistics.txt"
    k_len = None
    if stats_out.exists():
        try:
            with stats_out.open() as fh:
                for line in fh:
                    if line.startswith("k\t"):
                        try:
                            k_len = int(line.strip().split("\t", 1)[1])
                        except Exception:
                            k_len = None
                        break
        except Exception:
            k_len = None

    n_before = 0
    n_singleton_before = 0
    lens_pass = []
    lc_counts = Counter()
    removed_minlen = 0

    for _, seq in iter_fasta(in_fa):
        L = len(seq)
        n_before += 1
        if k_len is not None and L == k_len:
            n_singleton_before += 1
        if L < thr:
            removed_minlen += 1
            continue
        reason = lc_reason(seq)
        if reason is not None:
            lc_counts[reason] += 1
            continue
        lens_pass.append(L)

    removed_lc = sum(lc_counts.values())
    removed_total = removed_minlen + removed_lc
    n_after = len(lens_pass)
    total_bp_after = sum(lens_pass)
    n_singleton_after = sum(1 for L in lens_pass if k_len is not None and L == k_len)
    min_after = min(lens_pass) if lens_pass else 0
    max_after = max(lens_pass) if lens_pass else 0

    if removed_total == 0:
        pass
    else:
        tmp_out = res_dir / ".unitigs.filtered.tmp.fa"
        with open_text_write(tmp_out) as fout:
            idx = 0
            for _, seq in iter_fasta(in_fa):
                L = len(seq)
                if L < thr:
                    continue
                reason = lc_reason(seq)
                if reason is not None:
                    continue
                idx += 1
                fout.write(f">{idx} len_{L}\n")
                for i in range(0, L, 80):
                    fout.write(seq[i:i+80] + "\n")

        bak = None
        try:
            bak = res_dir / ".unitigs.fa.bak"
            if bak.exists():
                bak.unlink()
            in_fa.rename(bak)
        except Exception:
            bak = None
        shutil.move(str(tmp_out), str(res_dir / "unitigs.fa"))
        if bak is not None:
            try:
                bak.unlink()
            except Exception:
                pass

    mean_after = float(mean(lens_pass)) if lens_pass else 0.0
    median_after = float(median(lens_pass)) if lens_pass else 0.0
    n50_after, l50_after = n50_l50(lens_pass)

    try:
        cnt = Counter(lens_pass)
        lines = [
            "# Unitig statistics (unitigs.fa used for matching)",
            f"k\t{k_len if k_len is not None else 'NA'}",
            f"n_unitigs_incl_singletons\t{n_before}",
            f"unitigs_of_length_k\t{n_singleton_before if k_len is not None else 'NA'}",
            f"n_unitigs_for_matching\t{n_after}",
            f"n_unitigs\t{n_after}",
            f"total_bp\t{total_bp_after}",
            f"min_len\t{min_after if min_after is not None else 0}",
            f"max_len\t{max_after}",
            f"mean_len\t{mean_after}",
            f"median_len\t{median_after}",
            f"n50\t{n50_after}",
            f"l50\t{l50_after}",
            f"filter_min_length\t{thr}",
            f"filter_lc_homopolymer\toff",
            f"filter_lc_dinucleotide_repeat\toff",
            f"filter_lc_trinucleotide_repeat\toff",
            f"removed_min_length\t{removed_minlen}",
            f"removed_lc_total\t{removed_lc}",
            f"lc_reject_homopolymer\t{lc_counts.get('homopolymer', 0)}",
            f"lc_reject_dinucleotide_repeat\t{lc_counts.get('dinucleotide_repeat', 0)}",
            f"lc_reject_trinucleotide_repeat\t{lc_counts.get('trinucleotide_repeat', 0)}",
            f"removed_total\t{removed_total}",
            "note\tlength statistics & histogram reflect unitigs.fa used for matching",
            "",
            "# length_histogram (matching set)",
            "length\tcount",
        ]
        for L in sorted(cnt.keys()):
            lines.append(f"{L}\t{cnt[L]}")
        stats_out.write_text("\n".join(lines) + "\n")
    except Exception as exc:
        logger.warning(f"Could not rewrite unitig_statistics.txt: {exc}")

    summary = (
        f"Pass III summary: filtered_unitigs={n_after} "
        f"median_len_filtered={median_after:.2f} "
        f"mean_len_filtered={mean_after:.2f} "
        f"singletons_input={n_singleton_before if k_len is not None else 'NA'} "
        f"min_length={thr} removed_total={removed_total} "
        f"removed_minlen={removed_minlen} removed_lc={removed_lc} "
        f"(original_unitigs={n_before}; entropy filter disabled)"
    )
    try:
        with open(main_log, "a") as ml:
            ml.write(ts_line(summary) + "\n")
    except Exception as e:
        logger.warning(f"Could not append summary to main.log: {e}")

    logger.info(f"Filtered unitigs by min_length={thr}: before={n_before}, after={n_after}, removed_total={removed_total} (minlen={removed_minlen}, lc={removed_lc}, entropy=off)")

if __name__ == "__main__":
    sys.exit(main())
