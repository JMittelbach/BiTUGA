#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Create mapped and unmapped unitig ID lists from a unitig FASTA and a minimap2 PAF file. "
            "A unitig is considered mapped if its ID occurs at least once as qname in PAF column 1."
        )
    )
    ap.add_argument("--unitigs-fasta", type=Path, required=True, help="Input unitig FASTA file.")
    ap.add_argument("--paf", type=Path, required=True, help="Input minimap2 PAF file.")
    ap.add_argument("--mapped-ids-out", type=Path, required=True, help="Output TXT path for mapped IDs.")
    ap.add_argument("--unmapped-ids-out", type=Path, required=True, help="Output TXT path for unmapped IDs.")
    ap.add_argument(
        "--sort",
        choices=("numeric", "lex", "none"),
        default="numeric",
        help="Sort output IDs (default: numeric).",
    )
    ap.add_argument(
        "--strict",
        action="store_true",
        help="Fail if PAF contains query IDs that are not present in FASTA.",
    )
    return ap.parse_args()


def read_fasta_ids(path: Path) -> list[str]:
    ids: list[str] = []
    seen: set[str] = set()
    with path.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            uid = line[1:].strip().split()[0]
            if not uid:
                continue
            if uid in seen:
                raise ValueError(f"Duplicate FASTA ID: {uid}")
            seen.add(uid)
            ids.append(uid)
    return ids


def read_paf_query_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open() as fh:
        for lineno, line in enumerate(fh, start=1):
            raw = line.rstrip("\n")
            if not raw:
                continue
            cols = raw.split("\t")
            if len(cols) < 12:
                raise ValueError(f"Invalid PAF line {lineno}: expected >= 12 tab-separated columns")
            qid = cols[0]
            if qid:
                ids.add(qid)
    return ids


def key_numeric_or_lex(v: str):
    if v.isdigit():
        return (0, int(v))
    return (1, v)


def maybe_sort(ids: list[str], mode: str) -> list[str]:
    if mode == "none":
        return ids
    if mode == "lex":
        return sorted(ids)
    return sorted(ids, key=key_numeric_or_lex)


def write_ids(path: Path, ids: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for uid in ids:
            fh.write(uid + "\n")


def main() -> int:
    args = parse_args()

    if not args.unitigs_fasta.exists():
        print(f"[ERROR] FASTA not found: {args.unitigs_fasta}", file=sys.stderr)
        return 2
    if not args.paf.exists():
        print(f"[ERROR] PAF not found: {args.paf}", file=sys.stderr)
        return 2

    fasta_ids = read_fasta_ids(args.unitigs_fasta)
    fasta_set = set(fasta_ids)
    paf_ids = read_paf_query_ids(args.paf)

    unknown_in_paf = paf_ids - fasta_set
    if unknown_in_paf:
        msg = (
            f"[WARN] PAF has {len(unknown_in_paf)} query IDs not present in FASTA "
            f"(example: {next(iter(unknown_in_paf))})"
        )
        if args.strict:
            print(msg.replace("[WARN]", "[ERROR]"), file=sys.stderr)
            return 2
        print(msg, file=sys.stderr)

    mapped_set = fasta_set & paf_ids
    unmapped_set = fasta_set - mapped_set

    mapped_ids = maybe_sort(list(mapped_set), args.sort)
    unmapped_ids = maybe_sort(list(unmapped_set), args.sort)

    write_ids(args.mapped_ids_out, mapped_ids)
    write_ids(args.unmapped_ids_out, unmapped_ids)

    print(f"[INFO] FASTA IDs total      : {len(fasta_ids)}", file=sys.stderr)
    print(f"[INFO] PAF query IDs unique : {len(paf_ids)}", file=sys.stderr)
    print(f"[INFO] mapped IDs           : {len(mapped_ids)}", file=sys.stderr)
    print(f"[INFO] unmapped IDs         : {len(unmapped_ids)}", file=sys.stderr)
    print(f"[INFO] wrote mapped IDs to  : {args.mapped_ids_out}", file=sys.stderr)
    print(f"[INFO] wrote unmapped IDs to: {args.unmapped_ids_out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
