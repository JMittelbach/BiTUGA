#!/usr/bin/env python3
import sys, argparse, csv, os
from collections import Counter, defaultdict

def norm_trait(s: str) -> str:
    s = (s or "").strip().lower()
    s = " ".join(s.split())
    return s

def lev(a,b):
    la, lb = len(a), len(b)
    if la==0: return lb
    if lb==0: return la
    prev = list(range(lb+1))
    for i,ca in enumerate(a,1):
        cur = [i]
        for j,cb in enumerate(b,1):
            cost = 0 if ca==cb else 1
            cur.append(min(prev[j]+1, cur[j-1]+1, prev[j-1]+cost))
        prev = cur
    return prev[-1]

def detect_delim(path):
    with open(path, "r", newline="") as fh:
        head = fh.readline()
    if "\t" in head: return "\t"
    if ";" in head: return ";"
    if "," in head: return ","
    return None

ap = argparse.ArgumentParser(add_help=False)
ap.add_argument("--in", dest="inp", required=True)
ap.add_argument("--emit-sample-trait", dest="emit_trait")
ap.add_argument("--emit-sample-paths", dest="emit_paths")
ap.add_argument("--summary", action="store_true")
args = ap.parse_args()

meta_path = args.inp
out_trait = args.emit_trait
out_paths_dir = args.emit_paths
summary_only = bool(args.summary)

delim = detect_delim(meta_path)
rows = []
if delim:
    with open(meta_path, "r", newline="") as fh:
        rdr = csv.reader(fh, delimiter=delim)
        rows = [r for r in rdr]
else:
    with open(meta_path, "r", newline="") as fh:
        rows = [line.strip().split() for line in fh]

if not rows:
    print("ERROR: empty metadata", file=sys.stderr); sys.exit(1)

hdr = [c.strip() for c in rows[0]]
low = [c.lower() for c in hdr]
need = ["path","sample_id","trait"]
for k in need:
    if k not in low:
        print("ERROR: metadata requires headers: path, sample_id, trait", file=sys.stderr)
        sys.exit(1)
idx = {k: low.index(k) for k in need}

data = []
for r in rows[1:]:
    if not r: continue
    if len(r) < len(hdr): r = r + [""]*(len(hdr)-len(r))
    data.append([x for x in r])

traits_raw = [norm_trait(r[idx["trait"]]) for r in data if any(r)]
freq = Counter([t for t in traits_raw if t!=""])
if len(freq) < 2:
    print(f"ERROR: found only {len(freq)} unique trait value(s). Exactly two required.", file=sys.stderr)
    sys.exit(2)

seeds = [t for t,_ in freq.most_common(2)]
s1, s2 = seeds[0], seeds[1]

mapped_rows = []
bad = set()
for r in data:
    t = norm_trait(r[idx["trait"]])
    if t in (s1,s2):
        r[idx["trait"]] = t
        mapped_rows.append(r)
        continue
    d1 = lev(t, s1)
    d2 = lev(t, s2)
    if d1<=2 and d1<=d2:
        r[idx["trait"]] = s1
        mapped_rows.append(r)
    elif d2<=2:
        r[idx["trait"]] = s2
        mapped_rows.append(r)
    else:
        bad.add(t)

if bad:
    print("ERROR: Unmappable trait values detected. Only two consistent values are allowed (with typo tolerance).", file=sys.stderr)
    print("Unmappable values:", file=sys.stderr)
    for b in sorted(bad):
        print(f'  - "{b}"', file=sys.stderr)
    print(f"Detected target traits: {s1}, {s2}. Please correct {meta_path}.", file=sys.stderr)
    sys.exit(3)

sample_paths = defaultdict(list)
sample_trait = {}
for r in mapped_rows:
    p = r[idx["path"]].strip()
    sid = r[idx["sample_id"]].strip()
    tr = r[idx["trait"]].strip()
    if not p or not sid or not tr: continue
    sample_paths[sid].append(p)
    if sid not in sample_trait: sample_trait[sid] = tr

if summary_only:
    t_counts = Counter(sample_trait.values())
    if len(t_counts) != 2:
        print("ERROR: exactly two traits required after normalization.", file=sys.stderr)
        sys.exit(4)
    ordered = [t for t,_ in t_counts.most_common(2)]
    if t_counts[ordered[0]] == t_counts[ordered[1]] and ordered[0] > ordered[1]:
        ordered = [ordered[1], ordered[0]]
    t0, t1 = ordered[0], ordered[1]
    n0 = sum(1 for sid,tr in sample_trait.items() if tr == t0)
    n1 = sum(1 for sid,tr in sample_trait.items() if tr == t1)
    sys.stdout.write(f"{t0}\t{n0}\t{t1}\t{n1}\n")
    sys.exit(0)

with open(meta_path, "w", newline="") as fh:
    wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
    wr.writerow(hdr)
    wr.writerows(mapped_rows)

if out_trait:
    with open(out_trait, "w", newline="") as fh:
        for sid in sorted(sample_trait.keys()):
            fh.write(f"{sid}\t{sample_trait[sid]}\n")

if out_paths_dir:
    os.makedirs(out_paths_dir, exist_ok=True)
    for sid, paths in sample_paths.items():
        with open(os.path.join(out_paths_dir, sid), "w", newline="") as fh:
            for p in paths: fh.write(p+"\n")
