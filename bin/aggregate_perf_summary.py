#!/usr/bin/env python3
import argparse
import re
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple

PERF_RE = re.compile(
    r"^(?:\[\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2}\]\s+(?:\[[A-Z]+\]\s+)?)?"
    r"\[PERF\]\s+(?P<label>.+?):\s+runtime=(?P<runtime>\d{2}:\d{2}:\d{2}:\d{2}\.\d{3})\s+"
    r"max_rss=(?P<rss>[0-9.]+|NA)\s*(?:GB)?(?:\s+cpu_total_s=(?P<cpu>[0-9.NA]+)(?:\s+avg_cores=(?P<cores>[0-9.NA]+))?)?$"
)


def parse_runtime(runtime_str: str) -> float:
    """Return runtime in seconds from DD:HH:MM:SS.mmm."""
    dd, hh, mm, rest = runtime_str.split(":")
    ss, ms = rest.split(".")
    total_ms = (
        (int(dd) * 24 * 60 * 60 * 1000)
        + (int(hh) * 60 * 60 * 1000)
        + (int(mm) * 60 * 1000)
        + (int(ss) * 1000)
        + int(ms)
    )
    return total_ms / 1000.0


def parse_log(path: Path):
    entries = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        m = PERF_RE.match(line.strip())
        if not m:
            continue
        label = m.group("label")
        runtime_s = parse_runtime(m.group("runtime"))
        rss_raw = m.group("rss")
        rss_gb = float(rss_raw) if rss_raw != "NA" else None
        cpu_raw = m.group("cpu")
        cores_raw = m.group("cores")
        try:
            cpu_total = float(cpu_raw) if cpu_raw not in (None, "NA") else None
        except Exception:
            cpu_total = None
        try:
            avg_cores = float(cores_raw) if cores_raw not in (None, "NA") else None
        except Exception:
            avg_cores = None
        entries.append((label, runtime_s, rss_gb, cpu_total, avg_cores))
    return entries


def pass_label(label: str) -> str:
    """
    Normalize pass labels:
    - Map anything containing \"funnel plot\" or starting with \"Plotting\" to \"Plotting\"
    - Extract pass prefix (e.g., 'Pass IV') if present
    - Otherwise fall back to full label
    """
    low = label.lower()
    if "funnel plot" in low or low.startswith("plotting"):
        return "Plotting"
    m = re.search(r"(Pass\\s+[IVX]+)", label)
    return m.group(1) if m else label


def fmt_hms(seconds: float) -> str:
    s = int(seconds)
    ms = int(round((seconds - s) * 1000))
    m, s = divmod(s, 60)
    h, m = divmod(m, 60)
    d, h = divmod(h, 24)
    return f"{d:02d}:{h:02d}:{m:02d}:{s:02d}.{ms:03d}"


def main():
    ap = argparse.ArgumentParser(description="Aggregate PERF lines from BiTUGA main log.")
    ap.add_argument("--log", type=Path, default=Path("logging/main.log"),
                    help="Path to main log containing [PERF] lines (default: logging/main.log)")
    ap.add_argument("--out", type=Path, default=Path("logging/perf_summary.txt"),
                    help="Where to write summary (default: logging/perf_summary.txt)")
    args = ap.parse_args()

    if not args.log.is_file():
        raise SystemExit(f"Log not found: {args.log}")

    entries = parse_log(args.log)
    if not entries:
        raise SystemExit(f"No [PERF] lines found in {args.log}")

    by_pass = {}
    for label, runtime_s, rss_gb, cpu_total, _ in entries:
        key = pass_label(label)
        agg = by_pass.setdefault(key, {"runtime": 0.0, "rss": None, "cpu": 0.0, "cpu_present": False})
        if key == "Plotting":
            agg["runtime"] = runtime_s
            agg["cpu"] = cpu_total if cpu_total is not None else 0.0
            agg["cpu_present"] = cpu_total is not None
            if rss_gb is not None:
                agg["rss"] = rss_gb if agg["rss"] is None else max(agg["rss"], rss_gb)
        else:
            agg["runtime"] += runtime_s
            if rss_gb is not None:
                agg["rss"] = max(agg["rss"], rss_gb) if agg["rss"] is not None else rss_gb
            if cpu_total is not None:
                agg["cpu"] += cpu_total
                agg["cpu_present"] = True

    total_runtime = sum(vals["runtime"] for vals in by_pass.values())
    max_rss = max((vals["rss"] for vals in by_pass.values() if vals["rss"] is not None), default=None)
    total_cpu = sum((vals["cpu"] for vals in by_pass.values() if vals.get("cpu_present") and vals["cpu"] is not None and vals["cpu"] > 0), 0.0)
    cpu_present = any(vals.get("cpu_present") and vals["cpu"] is not None and vals["cpu"] > 0 for vals in by_pass.values())

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    prefix = lambda msg: f"[{ts}] [INFO] {msg}"
    lines = []
    lines.append(prefix("PERF SUMMARY"))
    lines.append(prefix(f"total_runtime={fmt_hms(total_runtime)}"))
    if max_rss is not None:
        lines.append(prefix(f"max_rss={max_rss:.3f} GB"))
    else:
        lines.append(prefix("max_rss=NA"))
    if cpu_present:
        lines.append(prefix(f"cpu_total_s={total_cpu:.3f}"))
    lines.append(prefix("by_pass:"))
    for key, vals in by_pass.items():
        rss_txt = f"{vals['rss']:.3f} GB" if vals["rss"] is not None else "NA"
        cpu_txt = ""
        if vals["cpu_present"] and vals["cpu"] is not None and vals["cpu"] > 0:
            cpu_txt = f", cpu_total_s={vals['cpu']:.3f}"
        lines.append(prefix(f"- {key}: runtime={fmt_hms(vals['runtime'])}, max_rss={rss_txt}{cpu_txt}"))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
