#!/usr/bin/env bash
set -euo pipefail

LABEL=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --label) LABEL="$2"; shift 2;;
    --) shift; break;;
    *) break;;
  esac
done
[[ -n "${LABEL}" ]] || { echo "[PERF] missing --label" >&2; exit 2; }
[[ $# -gt 0 ]] || { echo "[PERF] missing command" >&2; exit 2; }

ts(){ date '+%Y-%m-%d %H:%M:%S'; }

tmp="$(mktemp)"
start_ns="$(python3 - <<'PY'
import time; print(time.time_ns())
PY
)"

rc=0
if [[ "$(uname)" == "Darwin" ]]; then
  /usr/bin/time -l "$@" 2> "${tmp}" || rc=$?
else
  /usr/bin/time -v "$@" 2> "${tmp}" || rc=$?
fi
end_ns="$(python3 - <<'PY'
import time; print(time.time_ns())
PY
)"

wall_ns="$(( end_ns - start_ns ))"
wall_s="$(WALL_NS="$wall_ns" python3 - <<'PY'
import os
ns = int(os.environ["WALL_NS"])
print(ns / 1_000_000_000.0)
PY
)"

if [[ "$(uname)" == "Darwin" ]]; then
  rss_bytes="$(awk '/maximum resident set size|peak memory footprint/{val=$1} END{print val+0}' "${tmp}" 2>/dev/null || echo 0)"
  read usr_s sys_s <<<"$(python3 - "${tmp}" <<'PY'
import re, sys
txt = open(sys.argv[1], "r", errors="ignore").read()
def grab(pattern, default=0.0):
    m = re.search(pattern, txt)
    try:
        return float(m.group(1)) if m else default
    except Exception:
        return default
print(grab(r'([0-9.]+)\s+user'), grab(r'([0-9.]+)\s+sys'))
PY
)"
  rss_kb="$(( rss_bytes / 1024 ))"
  rss_gb_val="${rss_bytes}"
  rss_is_bytes="1"
else
  rss_kb="$(awk -F': *' '/Maximum resident set size/{print $2; exit}' "${tmp}" 2>/dev/null || echo 0)"
  read usr_s sys_s <<<"$(python3 - "${tmp}" <<'PY'
import re, sys
txt = open(sys.argv[1], "r", errors="ignore").read()
def grab(label):
    m = re.search(label + r':\s*([0-9.]+)', txt)
    try:
        return float(m.group(1)) if m else 0.0
    except Exception:
        return 0.0
print(grab("User time"), grab("System time"))
PY
)"
  rss_gb_val="${rss_kb}"
  rss_is_bytes="0"
fi

if [[ -z "${usr_s:-}" ]] || [[ -z "${sys_s:-}" ]]; then
  usr_s="nan"
  sys_s="nan"
fi
if [[ -z "${rss_kb:-}" ]] || [[ "${rss_kb}" == "0" && "${rss_is_bytes}" == "0" ]]; then
  rss_kb="nan"
fi

rm -f "${tmp}"

classic="$(python3 - "$wall_ns" "$usr_s" "$sys_s" "$rss_kb" "$rc" "$LABEL" <<'PY'
import sys
ns=int(sys.argv[1]); usr=float(sys.argv[2]); sys_t=float(sys.argv[3]); kb=int(float(sys.argv[4])); rc=int(sys.argv[5]); label=sys.argv[6]
runtime=ns//1_000_000_000
print(f"PERF {label}: runtime_s={runtime} user_s={usr} sys_s={sys_t} max_rss_kb={kb} rc={rc}")
PY
)"
pretty="$(python3 - "$wall_ns" "$rss_gb_val" "$rss_is_bytes" "$LABEL" "$usr_s" "$sys_s" "$wall_s" <<'PY'
import sys, math
ns=int(sys.argv[1]); val=sys.argv[2]; is_bytes=int(sys.argv[3]); label=sys.argv[4]; usr=float(sys.argv[5]); sys_t=float(sys.argv[6]); wall=float(sys.argv[7])
ms=ns//1_000_000
d=ms//86_400_000; ms%=86_400_000
h=ms//3_600_000;  ms%=3_600_000
m=ms//60_000;     ms%=60_000
s=ms//1000;       ms%=1000
try:
  x=float(val)
  gb = (x/1073741824.0) if is_bytes==1 else (x/1048576.0)
  mem = f"{gb:.3f} GB"
except:
  mem = "NA"
cpu_total = usr + sys_t if (not math.isnan(usr) and not math.isnan(sys_t)) else float('nan')
cpu_total_txt = f"{cpu_total:.3f}" if math.isfinite(cpu_total) else "NA"
print(f"[PERF] {label}: runtime={d:02d}:{h:02d}:{m:02d}:{s:02d}.{ms:03d}  max_rss={mem}  cpu_total_s={cpu_total_txt}")
PY
)"

printf "[%s] %s\n" "$(ts)" "${pretty}"
exit "${rc}"
