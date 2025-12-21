#!/usr/bin/env bash
set -euo pipefail

BIN_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
ROOT_DIR="$(cd -- "${BIN_DIR}/.." >/dev/null 2>&1 && pwd)"

usage() {
  cat >&2 <<'USAGE'
Usage:
  bin/build_groups_1.sh \
    [--metadata <metadata.tsv>] \
    [--k 31] [--ci 2] [--cx 4294967295] [--cs 2] \
    [--threads 8] [--mem-gb 8] \
    [--outdir <base_dir>] \
    [--kmc-bin <dir>] \
    [--gz-est-ratio <float>]   # Schätzwert für Entpack-Faktor bei *.gz (Default 3.5)
USAGE
  exit 1
}

METADATA="${ROOT_DIR}/metadata.tsv"
K=31; CI=2; CX=4294967295; CS=2
THREADS=8; MEM_GB=6
OUTDIR_BASE="${ROOT_DIR}"
KMC_BIN="${BIN_DIR}"
GZ_EST_RATIO="3.0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --metadata) METADATA="$2"; shift 2;;
    --k)        K="$2"; shift 2;;
    --ci)       CI="$2"; shift 2;;
    --cx)       CX="$2"; shift 2;;
    --cs)       CS="$2"; shift 2;;
    --threads)  THREADS="$2"; shift 2;;
    --mem-gb)   MEM_GB="$2"; shift 2;;
    --outdir)   OUTDIR_BASE="$2"; shift 2;;
    --kmc-bin)  KMC_BIN="$2"; shift 2;;
    --gz-est-ratio) GZ_EST_RATIO="$2"; shift 2;;
    -h|--help)  usage;;
    *) echo "[ERR] Unknown argument: $1" >&2; usage;;
  esac
done

case "${OUTDIR_BASE}" in /*) : ;; *) OUTDIR_BASE="${ROOT_DIR}/${OUTDIR_BASE}";; esac
case "${METADATA}"   in /*) : ;; *) METADATA="${ROOT_DIR}/${METADATA}";; esac
case "${KMC_BIN}"    in /*) : ;; *) KMC_BIN="${ROOT_DIR}/${KMC_BIN}";; esac

KMC="${KMC_BIN%/}/kmc"
KMCT="${KMC_BIN%/}/kmc_tools"
KMC_DUMP="${KMC_BIN%/}/kmc_dump"

RESULTS_DIR="${OUTDIR_BASE%/}/results"
POOLED_DIR="${RESULTS_DIR}/pooled_kmc_files"
STATS_DIR="${RESULTS_DIR}/statistics"
LOGROOT="${OUTDIR_BASE%/}/logging"
TMPDIR="${OUTDIR_BASE%/}/tmp"
LOGFILE="${LOGROOT}/build_pooled_kmc_files.log"

[[ -f "${METADATA}" ]] || { echo "[ERR] metadata not found: ${METADATA}" >&2; usage; }
mkdir -p "${RESULTS_DIR}" "${POOLED_DIR}" "${STATS_DIR}" "${LOGROOT}" "${TMPDIR}"

: > "${LOGFILE}"

ts(){ date '+%Y-%m-%d %H:%M:%S'; }
info(){ echo "[$(ts)] [INFO] $*"; }
warn(){ echo "[$(ts)] [WARN] $*" | tee -a "${LOGFILE}" >/dev/null; }
err(){  echo "[$(ts)] [ERROR] $*" 1>&2; echo "[$(ts)] [ERROR] $*" >> "${LOGFILE}"; }
trap 'err "Aborted at line ${LINENO} (Exit code $?)"' ERR

need_file(){ [[ -x "$1" ]] || { err "required binary not found or not executable: $1"; exit 127; }; }
need_file "${KMC}"
need_file "${KMCT}"

info "Starting KMC build & pooling"
info "Parameters: k=${K} ci=${CI} cx=${CX} cs=${CS} threads=${THREADS} mem_gb=${MEM_GB}"

SAMPLE_TRAIT_TSV="${TMPDIR}/sample_trait.tsv"
SAMPLE_PATHS_DIR="${TMPDIR}/sample_paths"
rm -f  "${SAMPLE_TRAIT_TSV}" 2>/dev/null || true
rm -rf "${SAMPLE_PATHS_DIR}" 2>/dev/null || true
mkdir -p "${SAMPLE_PATHS_DIR}"

NM="${BIN_DIR}/normalize_metadata_1.py"
[[ -x "${NM}" ]] || { err "normalize_metadata_1.py not executable: ${NM}"; exit 1; }

nm_out="$(
  set +e
  python3 "${NM}" \
    --in "${METADATA}" \
    --emit-sample-trait "${SAMPLE_TRAIT_TSV}" \
    --emit-sample-paths "${SAMPLE_PATHS_DIR}" \
    2>&1
  echo $? > "${TMPDIR}/.nm_rc"
  set -e
)"
nm_rc="$(cat "${TMPDIR}/.nm_rc" || echo 1)"
rm -f "${TMPDIR}/.nm_rc" 2>/dev/null || true

if [[ "${nm_rc}" -ne 0 ]]; then
  while IFS= read -r line; do
    [[ -n "${line}" ]] && err "${line}"
  done <<< "${nm_out}"
  exit "${nm_rc}"
fi

stat_bytes(){
  local f="$1"
  if [[ "$(uname)" == "Darwin" ]]; then
    stat -f%z -- "$f" 2>/dev/null || echo 0
  else
    stat -c%s -- "$f" 2>/dev/null || echo 0
  fi
}

bytes_human(){
  python3 - <<'PY' "$1"
import sys
b = float(sys.argv[1])
gb = b / (1024**3)
print(f"{gb:.2f} GB")
PY
}

total_files=0
total_bytes=0
gz_files=0
gz_bytes=0
est_uncompressed_bytes=0

while IFS= read -r pfile; do
  while IFS= read -r p; do
    [[ -z "$p" ]] && continue
    total_files=$(( total_files + 1 ))
    sz="$(stat_bytes "$p")"
    [[ "$sz" =~ ^[0-9]+$ ]] || sz=0
    total_bytes=$(( total_bytes + sz ))

    if [[ "$p" == *.gz ]]; then
      gz_files=$(( gz_files + 1 ))
      gz_bytes=$(( gz_bytes + sz ))
      est="$(python3 - <<PY "$sz" "$GZ_EST_RATIO"
import sys
sz=float(sys.argv[1]); ratio=float(sys.argv[2])
print(int(sz*ratio))
PY
)"
      est_uncompressed_bytes=$(( est_uncompressed_bytes + est ))
    else
      est_uncompressed_bytes=$(( est_uncompressed_bytes + sz ))
    fi
  done < "$pfile"
done < <(find "${SAMPLE_PATHS_DIR}" -type f -maxdepth 1 2>/dev/null | sort)

msg1="Input files: count=${total_files} (gz=${gz_files})"
msg2="Sizes: on-disk=$(bytes_human "${total_bytes}") ; gz-only=$(bytes_human "${gz_bytes}") ; est. uncompressed=$(bytes_human "${est_uncompressed_bytes}") [gz_est_ratio=${GZ_EST_RATIO}]"
info "${msg1}"
info "${msg2}"
{
  echo "[$(ts)] [INFO] ${msg1}"
  echo "[$(ts)] [INFO] ${msg2}"
} >> "${LOGFILE}"

TRAITS_LIST="${TMPDIR}/traits.list"
awk -F'\t' '{print $2}' "${SAMPLE_TRAIT_TSV}" | awk 'NF>0' | sort -u > "${TRAITS_LIST}"
trait_count="$(wc -l < "${TRAITS_LIST}" | tr -d '[:space:]')"
if [[ "${trait_count}" -ne 2 ]]; then err "metadata must contain exactly two traits; found: ${trait_count}"; exit 3; fi
trait1="$(sed -n '1p' "${TRAITS_LIST}")"
trait2="$(sed -n '2p' "${TRAITS_LIST}")"

build_sample_db() {
  sid="$1"; trait="$2"
  tmpd="${TMPDIR}/${sid}"
  mkdir -p "${tmpd}"
  : > "${tmpd}/.w"; rm -f "${tmpd}/.w" || { err "TMP not writable: ${tmpd}"; exit 2; }
  listfile="${tmpd}/${sid}.files.tsv"
  echo -e "path\tsample_id\ttrait" > "${listfile}"
  while IFS= read -r p; do
    [[ -n "$p" ]] && printf "%s\t%s\t%s\n" "$p" "$sid" "$trait"
  done < "${SAMPLE_PATHS_DIR}/${sid}" >> "${listfile}"

  echo "[create] ${BIN_DIR}/create_kmc3_sample_1.sh --sample ${sid}" >> "${LOGFILE}"
  info "Create sample ${sid}"
  "${BIN_DIR}/create_kmc3_sample_1.sh" \
      --files-tsv "${listfile}" \
      --sample    "${sid}" \
      --k         "${K}" \
      --ci        "${CI}" \
      --cx        "${CX}" \
      --cs        "${CS}" \
      --threads   "${THREADS}" \
      --mem-gb    "${MEM_GB}" \
      --tmpdir    "${tmpd}" \
      --logfile   "${LOGFILE}" \
      --kmc-bin   "${KMC_BIN}"
  echo "[create] DONE ${sid}" >> "${LOGFILE}"
  info "Create done ${sid}"

  rm -f "${listfile}"
}

merge_and_delete() {
  sid="$1"; trait="$2"
  db_dir="${TMPDIR}/${sid}"
  db="${db_dir}/${sid}"
  if [[ ! -f "${db}.kmc_pre" || ! -f "${db}.kmc_suf" ]]; then
    warn "DB missing for ${sid} in TMP — skipping merge"
    return
  fi
  info "Merging into pool: ${sid} ${trait}"
  "${BIN_DIR}/update_pool_1.sh" --sample-db "${db}" --trait "${trait}" --outdir "${RESULTS_DIR}" --kmc-bin "${KMC_BIN}" >> "${LOGFILE}" 2>&1
  rm -f "${db}.kmc_pre" "${db}.kmc_suf" || true
  rmdir "${db_dir}" 2>/dev/null || true
}

for trait in "${trait1}" "${trait2}"; do
  info "Processing trait=${trait}"
  while IFS=$'\t' read -r sid t; do
    [[ -z "${sid:-}" ]] && continue
    [[ "${t:-}" == "${trait}" ]] || continue
    build_sample_db "${sid}" "${t}"
    merge_and_delete "${sid}" "${t}"
  done < "${SAMPLE_TRAIT_TSV}"
  info "Pool completed: ${trait} → pooled_${trait}.kmc_pre / .kmc_suf"
done

info "Aggregating stats"
"${BIN_DIR}/aggregate_kmc_stats_and_totals_1.sh" --metadata "${METADATA}" --outdir "${RESULTS_DIR}" --logfile "${LOGFILE}" >> "${LOGFILE}" 2>&1
info "Done"

if [[ -d "${TMPDIR}" ]]; then
  rm -rf "${TMPDIR}" || warn "Failed to delete TMPDIR"
fi
