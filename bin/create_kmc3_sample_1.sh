#!/usr/bin/env bash
set -euo pipefail
usage(){ echo "Usage: $0 --files-tsv <tsv> --sample <id> --k <int> --ci <int> --cx <int> --cs <int> --threads <int> --mem-gb <int> --tmpdir <dir> --logfile <path> --kmc-bin <dir>"; exit 1; }

FILES_TSV="" ; SAMPLE="" ; K=31
CI=1 ; CX=4294967295 ; CS=4294967295
THREADS=4 ; MEM_GB=4
TMPDIR="./tmp"
LOGFILE=""
KMC_BIN=""

BIN_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
ROOT_DIR="$(cd -- "${BIN_DIR}/.." >/dev/null 2>&1 && pwd)"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --files-tsv) FILES_TSV="$2"; shift 2;;
    --sample)    SAMPLE="$2"; shift 2;;
    --k)         K="$2"; shift 2;;
    --ci)        CI="$2"; shift 2;;
    --cx)        CX="$2"; shift 2;;
    --cs)        CS="$2"; shift 2;;
    --threads)   THREADS="$2"; shift 2;;
    --mem-gb)    MEM_GB="$2"; shift 2;;
    --tmpdir)    TMPDIR="$2"; shift 2;;
    --logfile)   LOGFILE="$2"; shift 2;;
    --kmc-bin)   KMC_BIN="$2"; shift 2;;
    -h|--help)   usage;;
    *) echo "Unknown arg: $1" >&2; usage;;
  esac
done

[[ -f "$FILES_TSV" && -n "$SAMPLE" ]] || usage
[[ -n "$LOGFILE" ]] || LOGFILE="${ROOT_DIR}/logging/build_pooled_kmc_files.log"
mkdir -p "$(dirname "${LOGFILE}")"

[[ -n "$KMC_BIN" ]] || KMC_BIN="${BIN_DIR}"
case "${KMC_BIN}" in /*) : ;; *) KMC_BIN="${ROOT_DIR}/${KMC_BIN}";; esac

KMC="${KMC_BIN%/}/kmc"
[[ -x "${KMC}" ]] || { echo "[create_kmc3_sample] ERR: kmc not executable at ${KMC}" >> "${LOGFILE}"; exit 127; }

mkdir -p "${TMPDIR}"

DB_RAW="${TMPDIR%/}/${SAMPLE}_raw"
DB_FINAL="${TMPDIR%/}/${SAMPLE}"
LIST="${TMPDIR%/}/${SAMPLE}.files.txt"

awk -F'\t' -v s="${SAMPLE}" '
NR==1{ for(i=1;i<=NF;i++) h[$i]=i; next }
$h["sample_id"]==s { print $h["path"] }
' "${FILES_TSV}" > "${LIST}" || true
[[ -s "${LIST}" ]] || { echo "[create_kmc3_sample] No files for sample=${SAMPLE}" >> "${LOGFILE}"; exit 1; }

missing=0
while IFS= read -r p; do
  [[ -z "$p" ]] && continue
  [[ -f "$p" ]] || { echo "[create_kmc3_sample] Missing input: $p" >> "${LOGFILE}"; missing=$((missing+1)); }
done < "${LIST}"
[[ "$missing" -eq 0 ]] || { echo "[create_kmc3_sample] Abort: ${missing} missing files for ${SAMPLE}" >> "${LOGFILE}"; exit 2; }

first=$(head -n1 "${LIST}")
mode="-fq"
case "${first##*.}" in
  fa|fasta|fa.gz|fasta.gz) mode="-fa";;
  fq|fastq|fq.gz|fastq.gz) mode="-fq";;
  gz) base=$(basename "$first" .gz); case "${base##*.}" in fa|fasta) mode="-fa";; *) mode="-fq";; esac;;
  *) mode="-fq";;
esac

STDBUF=""
if command -v stdbuf >/dev/null 2>&1; then
  STDBUF="stdbuf -oL -eL"
fi

PROG_TMP="$(mktemp "${TMPDIR%/}/kmc_progress.${SAMPLE}.XXXXXX")"
cleanup(){ rm -f "${PROG_TMP}" 2>/dev/null || true; }
trap cleanup EXIT

printf "%s\n" "**********************" >> "${LOGFILE}"

set +e
${STDBUF} "${KMC}" -k"${K}" -m"${MEM_GB}" -t"${THREADS}" -ci"${CI}" -cx"${CX}" -cs"${CS}" \
    ${mode} @"${LIST}" "${DB_RAW}" "${TMPDIR}" 1>>"${LOGFILE}" 2>>"${PROG_TMP}"
rc=$?
set -e

if [[ -s "${PROG_TMP}" ]]; then
  tr '\r' '\n' < "${PROG_TMP}" | sed -E '/^[[:space:]]*$/d' >> "${LOGFILE}"
fi

if [[ $rc -ne 0 ]]; then
  echo "[create_kmc3_sample] kmc failed (rc=${rc})" >> "${LOGFILE}"
  exit $rc
fi

mv -f "${DB_RAW}.kmc_pre" "${DB_FINAL}.kmc_pre"
mv -f "${DB_RAW}.kmc_suf" "${DB_FINAL}.kmc_suf"
rm -f "${LIST}"

echo "[create_kmc3_sample] DONE: ${DB_FINAL}.kmc_pre / .kmc_suf" >> "${LOGFILE}"
