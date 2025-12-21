#!/usr/bin/env python3
import argparse, gzip, logging, shutil, subprocess, sys, tempfile, time
from pathlib import Path
from statistics import mean, median
from collections import Counter

def ensure_dir(p: Path): p.mkdir(parents=True, exist_ok=True)
def open_text(p: Path): return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "rt")

def iter_fasta(p: Path):
    with open_text(p) as fh:
        h=None; buf=[]
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(buf)
                h = line[1:].strip(); buf=[]
            else:
                buf.append(line.strip())
        if h is not None:
            yield h, "".join(buf)

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

def compute_stats(fa: Path, k: int):
    lens = []
    n_single = 0
    for _, seq in iter_fasta(fa):
        L = len(seq)
        lens.append(L)
        if L == k:
            n_single += 1
    n = len(lens)
    if n == 0:
        return {
            "k": k,
            "n_incl_singletons": 0,
            "unitigs_of_length_k": 0,
            "n_unitigs_for_matching": 0,
            "n_unitigs": 0,
            "total_bp": 0,
            "min_len": 0,
            "max_len": 0,
            "mean_len": 0.0,
            "median_len": 0.0,
            "n50": 0,
            "l50": 0,
            "lengths": [],
        }
    total = sum(lens)
    n50, l50 = n50_l50(lens)
    return {
        "k": k,
        "n_incl_singletons": n,
        "unitigs_of_length_k": n_single,
        "n_unitigs_for_matching": n,
        "n_unitigs": n,
        "total_bp": total,
        "min_len": min(lens),
        "max_len": max(lens),
        "mean_len": float(mean(lens)),
        "median_len": float(median(lens)),
        "n50": n50,
        "l50": l50,
        "lengths": lens,
    }

def write_stats(stats: dict, outp: Path):
    cnt = Counter(stats["lengths"])
    with open(outp, "w") as fh:
        fh.write("# Unitig statistics (unitigs.fa used for matching)\n")
        fh.write(f"k\t{stats['k']}\n")
        fh.write(f"n_unitigs_incl_singletons\t{stats['n_incl_singletons']}\n")
        fh.write(f"unitigs_of_length_k\t{stats['unitigs_of_length_k']}\n")
        fh.write(f"n_singleton_k\t{stats['unitigs_of_length_k']}\n")
        fh.write(f"n_unitigs_for_matching\t{stats['n_unitigs_for_matching']}\n")
        fh.write(f"n_unitigs\t{stats['n_unitigs_for_matching']}\n")
        fh.write(f"total_bp\t{stats['total_bp']}\n")
        fh.write(f"min_len\t{stats['min_len']}\n")
        fh.write(f"max_len\t{stats['max_len']}\n")
        fh.write(f"mean_len\t{stats['mean_len']}\n")
        fh.write(f"median_len\t{stats['median_len']}\n")
        fh.write(f"n50\t{stats['n50']}\n")
        fh.write(f"l50\t{stats['l50']}\n")
        fh.write("note\tlength statistics & histogram reflect unitigs.fa used for matching\n")
        fh.write("\n# length_histogram (matching set)\nlength\tcount\n")
        for L in sorted(cnt.keys()):
            fh.write(f"{L}\t{cnt[L]}\n")

def setup_logger(log_dir: Path):
    ensure_dir(log_dir)
    logger = logging.getLogger("build_unitigs")
    if logger.handlers:
        logger.handlers.clear()
    logger.setLevel(logging.INFO)
    logger.propagate = False
    fmt=logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh=logging.FileHandler(log_dir/"build_unitigs.log", mode="a")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    return logger

def ts_line(msg: str) -> str:
    return f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] [INFO] {msg}"

def main():
    ap=argparse.ArgumentParser(description=("Build unitigs from candidate_kmers.fasta(.gz) using BCALM (default k=31). Writes only results/unitigs.fa + one statistics file; all BCALM intermediates go to a temp workdir and are removed."))
    ap.add_argument("--base-dir", default=None)
    ap.add_argument("--candidate","-c", default=None)
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--abundance-min", type=int, default=1)
    ap.add_argument("--minimizer-size", type=int, default=None)
    ap.add_argument("--max-memory-gb", type=int, default=None)
    args=ap.parse_args()

    script_path=Path(__file__).resolve()
    bin_dir=script_path.parent
    repo_root=bin_dir.parent
    base_dir=Path(args.base_dir).resolve() if args.base_dir else repo_root

    results_dir=base_dir/"results"
    logging_dir=base_dir/"logging"
    stats_dir=results_dir/"statistics"
    ensure_dir(results_dir); ensure_dir(logging_dir); ensure_dir(stats_dir)

    logger=setup_logger(logging_dir)
    main_log_path = logging_dir/"main.log"

    bcalm_path=bin_dir/"bcalm"
    if not bcalm_path.exists():
        print(f"Error: {bcalm_path} not found.", file=sys.stderr); sys.exit(2)

    cand = Path(args.candidate).resolve() if args.candidate else (results_dir/"candidate_kmers.fasta")
    if not cand.exists():
        alt = results_dir/"candidate_kmers.fasta.gz"
        if alt.exists():
            cand = alt
        else:
            logger.error(f"Input not found: {cand} (or {results_dir/'candidate_kmers.fasta.gz'})")
            sys.exit(3)

    final_unitigs = results_dir/"unitigs.fa"
    bcalm_log = logging_dir/"build_unitigs.log"

    temp_workdir = Path(tempfile.mkdtemp(prefix="unitigs_build_", dir=str(results_dir)))
    out_prefix = temp_workdir/"unitigs"

    cmd=[str(bcalm_path),
         "-in", str(cand),
         "-kmer-size", str(args.k),
         "-abundance-min", str(args.abundance_min),
         "-nb-cores", str(args.threads),
         "-out", str(out_prefix)]
    if args.minimizer_size is not None:
        cmd+=["-minimizer-size", str(args.minimizer_size)]
    if args.max_memory_gb is not None:
        cmd+=["-max-memory", str(args.max_memory_gb*1024)]

    logger.info(f"Using temp workdir: {temp_workdir}")
    logger.info("CMD: "+" ".join(cmd))
    try:
        with open(bcalm_log,"a") as so:
            so.write("CMD: "+" ".join(cmd)+"\n"); so.flush()
            proc=subprocess.run(cmd, stdout=so, stderr=subprocess.DEVNULL)
        if proc.returncode!=0:
            logger.error(f"BCALM exited with code {proc.returncode}. See {bcalm_log}")
            sys.exit(proc.returncode)

        produced=None
        for suf in (".unitigs.fa",".unitigs.fa.gz",".unitigs.fasta",".unitigs.fasta.gz"):
            p=Path(str(out_prefix)+suf)
            if p.exists():
                produced=p; break
        if produced is None:
            logger.error("BCALM output not found in temp workdir.")
            sys.exit(4)

        if produced.suffix==".gz":
            with gzip.open(produced,"rb") as fin, open(final_unitigs,"wb") as fout:
                shutil.copyfileobj(fin,fout)
        else:
            shutil.copyfile(produced, final_unitigs)
        logger.info(f"Final unitigs written to {final_unitigs}")

        stats = compute_stats(final_unitigs, args.k)
        stats_txt = stats_dir/"unitig_statistics.txt"
        write_stats(stats, stats_txt)
        logger.info(f"Statistics written: {stats_txt}")

        summary = (
            f"Pass III summary: unitigs_matching={stats['n_unitigs_for_matching']}  "
            f"median_len={stats['median_len']:.2f}  "
            f"singletons_in_input(k={args.k})={stats['unitigs_of_length_k']}"
        )
        try:
            with open(main_log_path, "a") as ml:
                ml.write(ts_line(summary) + "\n")
        except Exception as e:
            logger.warning(f"Could not append summary to main.log: {e}")

    finally:
        try:
            shutil.rmtree(temp_workdir)
            logger.info(f"Cleaned temp workdir: {temp_workdir}")
        except Exception as e:
            logger.warning(f"Could not remove temp workdir {temp_workdir}: {e}")

    logger.info(f"Done. (bcalm log: {bcalm_log})")

if __name__=="__main__":
    main()
