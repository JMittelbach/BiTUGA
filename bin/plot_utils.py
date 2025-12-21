#!/usr/bin/env python3
from pathlib import Path
from typing import Optional


def read_k(run_root: Path) -> Optional[int]:
    """Try to read k-mer length k from various stats files."""
    candidates = [
        run_root / "results/statistics/unitig_statistics.txt",
        run_root / "results/statistics/kmer_candidates.txt",
        run_root / "results/statistics/unitig_association_stats.txt",
        run_root / "results/statistics/unitig_matches.txt",
    ]
    for path in candidates:
        if not path.is_file():
            continue
        try:
            with path.open() as fh:
                for line in fh:
                    if line.startswith("k\t"):
                        parts = line.strip().split("\t")
                        if len(parts) >= 2:
                            return int(float(parts[1]))
                        try:
                            return int(float(line.strip().split()[-1]))
                        except Exception:
                            continue
        except Exception:
            continue
    return None
