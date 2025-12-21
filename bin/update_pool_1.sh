#!/usr/bin/env bash
set -euo pipefail
usage(){ echo "Usage: $0 --sample-db <path_no_ext> --trait <name> --outdir <results_dir> --kmc-bin <dir>"; exit 1; }

SAMPLE_DB="" ; TRAIT="" ; OUTDIR="results" ; KMC_BIN=""
while [[ $# -gt 0 ]]; do case "$1" in
  --sample-db) SAMPLE_DB="$2"; shift 2;;
  --trait)     TRAIT="$2"; shift 2;;
  --outdir)    OUTDIR="$2"; shift 2;;
  --kmc-bin)   KMC_BIN="$2"; shift 2;;
  -h|--help)   usage;;
  *) echo "Unknown arg: $1"; usage;;
esac; done

[[ -n "$SAMPLE_DB" && -n "$TRAIT" ]] || usage
[[ -n "${KMC_BIN}" ]] || KMC_BIN="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
case "${KMC_BIN}" in /*) : ;; *) KMC_BIN="$(cd -- "$(dirname "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)/${KMC_BIN}";; esac
KMCT="${KMC_BIN%/}/kmc_tools"
[[ -x "${KMCT}" ]] || { echo "[update_pool] ERR: kmc_tools not executable at ${KMCT}"; exit 127; }

POOLED_DIR="${OUTDIR%/}/pooled_kmc_files"
mkdir -p "${POOLED_DIR}"
POOL="${POOLED_DIR}/pooled_${TRAIT}"

if [[ ! -f "${POOL}.kmc_pre" || ! -f "${POOL}.kmc_suf" ]]; then
  cp "${SAMPLE_DB}.kmc_pre" "${POOL}.kmc_pre"
  cp "${SAMPLE_DB}.kmc_suf" "${POOL}.kmc_suf"
  exit 0
fi

TMP="${POOL}_tmp"
"${KMCT}" simple "${POOL}" "${SAMPLE_DB}" union "${TMP}"
mv "${TMP}.kmc_pre" "${POOL}.kmc_pre"
mv "${TMP}.kmc_suf" "${POOL}.kmc_suf"
