#!/usr/bin/env python3
import argparse
import sys
import re
from pathlib import Path
from collections import Counter, defaultdict

UNCOMPRESSED_FACTOR = 3.0

def parse_args():
    script_path = Path(__file__).resolve()
    root_dir = script_path.parent.parent
    default_output = root_dir / "metadata.tsv"
    default_traits = root_dir / "trait_info.tsv"

    p = argparse.ArgumentParser(
        description="Generate metadata and size summary for FASTQ files."
    )
    p.add_argument(
        "-o", "--output",
        default=str(default_output),
        help=f"Path to metadata TSV output file (default: {default_output})"
    )
    p.add_argument(
        "-t", "--traits",
        default=str(default_traits),
        help=f"Path to trait info file (TSV/CSV/TXT). Default: {default_traits}"
    )
    p.add_argument(
        "--exclude",
        nargs="*",
        default=[],
        help="Sample IDs to exclude completely (space- or comma-separated list)."
    )
    p.add_argument(
        "-d", "--dirs",
        nargs="+",
        action="append",
        required=False,
        help="One or more directories containing FASTQ/FASTQ.GZ files. (Deprecated, use --input-dir.)"
    )
    p.add_argument(
        "--input-dir",
        nargs="+",
        action="append",
        required=False,
        help="One or more directories containing FASTQ/FASTQ.GZ files (can be given multiple times)."
    )
    args = p.parse_args()

    dirs = []
    for group in args.dirs or []:
        dirs.extend(group)
    for group in args.input_dir or []:
        dirs.extend(group)
    if not dirs:
        p.error("At least one --input-dir (or --dirs) is required.")
    args.dirs = dirs
    return args

def iter_fastq_files(base_dirs):
    seen = set()
    patterns = ("*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz")
    for d in base_dirs:
        dpath = Path(d)
        if not dpath.is_dir():
            print(f"Warning: {dpath} is not a directory, skipping", file=sys.stderr)
            continue
        for pattern in patterns:
            for f in dpath.glob(pattern):
                rf = f.resolve()
                if rf in seen:
                    continue
                seen.add(rf)
                yield rf

def guess_sample_id(filename: str, known_ids) -> str:
    """
    Try to infer sample_id from filename using known IDs from trait file.
    Prefer exact substring matches (case-insensitive) over fallback split-on-underscore.
    """
    base = Path(filename).name
    stem = re.sub(r"\.(fastq|fq)(\.gz)?$", "", base, flags=re.IGNORECASE)
    stem_lower = stem.lower()
    matches = []
    for kid in known_ids:
        kid_lower = kid.lower()
        if kid_lower in stem_lower:
            matches.append(kid)
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        chosen = max(matches, key=len)
        print(f"Warning: multiple sample_id matches in {base}: {matches} ; choosing {chosen}", file=sys.stderr)
        return chosen
    return base.split("_", 1)[0]

def sample_sort_key(sample_id: str):
    m = re.match(r"^[A-Za-z]*([0-9]+)$", sample_id)
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            return sample_id
    return sample_id

def clean_trait(raw: str) -> str:
    s = raw.strip().lower()
    s = re.sub(r"[^a-z0-9äöüß]", "", s)
    return s

def levenshtein(a: str, b: str) -> int:
    if a == b:
        return 0
    la, lb = len(a), len(b)
    if la == 0:
        return lb
    if lb == 0:
        return la
    prev_row = list(range(lb + 1))
    for i in range(1, la + 1):
        cur_row = [i] + [0] * lb
        for j in range(1, lb + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            cur_row[j] = min(
                prev_row[j] + 1,
                cur_row[j - 1] + 1,
                prev_row[j - 1] + cost,
            )
        prev_row = cur_row
    return prev_row[-1]

def normalize_sex_label(raw: str) -> str | None:
    t = raw.strip().lower()
    t_letters = re.sub(r"[^a-zäöüß]", "", t)
    male_vals = {"m", "male", "mann", "maennlich", "männlich"}
    female_vals = {"f", "female", "frau", "weiblich"}
    if t_letters in male_vals:
        return "male"
    if t_letters in female_vals:
        return "female"
    if levenshtein(t_letters, "male") <= 1 or levenshtein(t_letters, "mann") <= 1:
        return "male"
    if levenshtein(t_letters, "female") <= 1 or levenshtein(t_letters, "weiblich") <= 2:
        return "female"
    return None

def build_trait_mapping(pairs):
    raw_traits = [tr for _, tr in pairs if tr.strip()]
    if not raw_traits:
        return {}

    base_to_raw_counts = defaultdict(Counter)
    for _, tr in pairs:
        base = clean_trait(tr)
        if not base:
            base = tr.strip().lower()
        base_to_raw_counts[base][tr] += 1

    def base_sex(base: str) -> str | None:
        for raw in base_to_raw_counts[base].keys():
            return normalize_sex_label(raw)
        return None

    bases = list(base_to_raw_counts.keys())
    clusters = []

    for base in bases:
        placed = False
        sex_base = base_sex(base)
        for cluster in clusters:
            rep = cluster[0]
            sex_rep = base_sex(rep)

            if sex_base is not None and sex_rep is not None and sex_base != sex_rep:
                continue

            if base == rep or levenshtein(base, rep) <= 2 or sorted(base) == sorted(rep):
                cluster.append(base)
                placed = True
                break
        if not placed:
            clusters.append([base])

    raw_to_canonical = {}
    canonical_labels = []

    for cluster in clusters:
        combined = Counter()
        for base in cluster:
            combined.update(base_to_raw_counts[base])
        canonical_raw, _ = combined.most_common(1)[0]
        sex_norm = normalize_sex_label(canonical_raw)
        if sex_norm is not None:
            canonical = sex_norm
        else:
            canonical = canonical_raw.strip()
        canonical_labels.append(canonical)
        for base in cluster:
            for raw in base_to_raw_counts[base].keys():
                raw_to_canonical[raw] = canonical

    unique_labels = sorted(set(canonical_labels))
    print(f"Detected {len(unique_labels)} trait value(s): {', '.join(unique_labels)}", file=sys.stderr)
    if len(unique_labels) > 2:
        print("Warning: more than two distinct trait values detected.", file=sys.stderr)
    return raw_to_canonical

def read_traits(path: Path):
    if not path.is_file():
        print(f"Error: trait file {path} not found", file=sys.stderr)
        sys.exit(1)
    pairs = []
    splitter = re.compile(r"[,\t;'\-]+|\s+")
    header_seen = False
    sid_idx = 0
    tr_idx = 1
    allowed_headers = {"id", "sample", "sample_id", "sampleid"}
    with path.open("r", encoding="utf-8") as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in splitter.split(line) if p.strip()]
            if not parts:
                continue
            if not header_seen:
                cleaned = [re.sub(r"[^a-z0-9]", "", p.lower()) for p in parts]
                best_idx = None
                best_dist = 3
                for idx, c in enumerate(cleaned):
                    for target in allowed_headers:
                        dist = levenshtein(c, target)
                        if dist < best_dist and dist <= 2:
                            best_dist = dist
                            best_idx = idx
                if best_idx is not None:
                    sid_idx = best_idx
                    tr_idx = 0 if sid_idx != 0 else (1 if len(parts) > 1 else 0)
                    header_seen = True
                    continue
                header_seen = True
            if len(parts) < 2:
                continue
            if sid_idx >= len(parts) or tr_idx >= len(parts):
                continue
            sid = parts[sid_idx].strip()
            trait_raw = parts[tr_idx].strip()
            if not sid:
                continue
            pairs.append((sid, trait_raw))
    raw_to_canonical = build_trait_mapping(pairs)
    traits = {}
    for sid, trait_raw in pairs:
        traits[sid] = raw_to_canonical.get(trait_raw, trait_raw.strip())
    print(f"Loaded {len(traits)} trait entries from {Path(path).name}", file=sys.stderr)
    return traits

def parse_excludes(raw_list):
    ids = set()
    for token in raw_list:
        for part in token.replace(",", " ").split():
            part = part.strip()
            if part:
                ids.add(part)
    return ids

def main():
    args = parse_args()
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    trait_path = Path(args.traits)
    traits = read_traits(trait_path)
    known_ids = set(traits.keys())

    excluded_ids = parse_excludes(args.exclude)

    rows = []
    sizes = {}
    missing_trait_ids = set()

    fastq_files = list(iter_fastq_files(args.dirs))
    total_bytes = sum(f.stat().st_size for f in fastq_files)
    total_gb = total_bytes / 1e9
    print(f"[INFO] Scanning {len(fastq_files)} FASTQ files (~{total_gb:.2f} GB compressed) to "
          "collect paths and sizes...", file=sys.stderr)

    for f in fastq_files:
        sample_id = guess_sample_id(f.name, known_ids)
        if sample_id in excluded_ids:
            continue
        trait = traits.get(sample_id, "NA")
        if trait == "NA":
            missing_trait_ids.add(sample_id)

        size_bytes = f.stat().st_size
        is_gz = "".join(f.suffixes).endswith(".gz")

        if sample_id not in sizes:
            sizes[sample_id] = {
                "n_files": 0,
                "compressed_bytes": 0.0,
                "uncompressed_est_bytes": 0.0,
            }

        sizes[sample_id]["n_files"] += 1
        sizes[sample_id]["compressed_bytes"] += size_bytes
        if is_gz:
            sizes[sample_id]["uncompressed_est_bytes"] += size_bytes * UNCOMPRESSED_FACTOR
        else:
            sizes[sample_id]["uncompressed_est_bytes"] += size_bytes

        rows.append((str(f), sample_id, trait))

    rows.sort(key=lambda r: (sample_sort_key(r[1]), r[0]))

    with out_path.open("w", encoding="utf-8") as out_f:
        out_f.write("path\tsample_id\ttrait\n")
        for path_str, sample_id, trait in rows:
            out_f.write(f"{path_str}\t{sample_id}\t{trait}\n")

    sizes_path = out_path.with_name(out_path.stem + "_sizes.tsv")
    with sizes_path.open("w", encoding="utf-8") as s_f:
        s_f.write("sample_id\tn_files\ttotal_GB\test_uncompressed_GB\n")
        total_files = 0
        total_comp_bytes = 0.0
        total_uncomp_bytes = 0.0

        for sample_id in sorted(sizes.keys(), key=sample_sort_key):
            info = sizes[sample_id]
            n_files = info["n_files"]
            comp_gb = info["compressed_bytes"] / 1e9
            uncomp_gb = info["uncompressed_est_bytes"] / 1e9
            total_files += n_files
            total_comp_bytes += info["compressed_bytes"]
            total_uncomp_bytes += info["uncompressed_est_bytes"]
            s_f.write(f"{sample_id}\t{n_files}\t{comp_gb:.3f}\t{uncomp_gb:.3f}\n")

        total_comp_gb = total_comp_bytes / 1e9
        total_uncomp_gb = total_uncomp_bytes / 1e9
        s_f.write(f"TOTAL\t{total_files}\t{total_comp_gb:.3f}\t{total_uncomp_gb:.3f}\n")

    print(f"Metadata written to: {out_path.name} ({len(rows)} rows)", file=sys.stderr)
    print(f"Size summary written to: {sizes_path.name} ({len(sizes)} samples)", file=sys.stderr)
    if missing_trait_ids:
        miss_str = ", ".join(sorted(missing_trait_ids, key=sample_sort_key))
        print(f"Warning: missing traits for sample IDs: {miss_str}", file=sys.stderr)

if __name__ == "__main__":
    main()
