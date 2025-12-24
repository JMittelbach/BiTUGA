#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
PASS0_BIN="${ROOT_DIR}/bin/make_metadata.py"
BUILD="${ROOT_DIR}/bin/build_groups_1.sh"
PASS2_BIN="${ROOT_DIR}/bin/merge2stats"
NM="${ROOT_DIR}/bin/normalize_metadata_1.py"
PERF="${ROOT_DIR}/bin/perf_run.sh"
PASS3_PY="${ROOT_DIR}/bin/build_unitigs_3.py"
PASS4_BIN="${ROOT_DIR}/bin/orchestrator_4.py"
PASS5_BIN="${ROOT_DIR}/bin/orchestrator_5.py"
QSPLIT_DEFAULT="262144"

usage() {
  cat >&2 <<'USAGE'
Usage: ./BiTUGA.sh [options]

Key options:
  Option                 Default         Description
  --input-dir <dir>      -               One or more FASTQ directories (repeat the flag for multiple dirs)
  --trait-info <file>    trait_info.tsv  Trait info table (tsv/csv/txt; header with ID and trait column; delim: tab/comma/space)
  --out-dir <dir>        current dir     Output and working directory
  --threads <int>        12              Max CPU threads to use
  --mem-gb <int>         16              Max RAM (GB) to use
  --k <int>              31              k-mer length
  --prev-min <float>     0.25            Keep k-mers present in at least this fraction of samples
  --prev-max <float>     0.75            Drop k-mers present in more than this fraction of samples
  --group-contrast       off             Shortcut: set group-min=0.50 and group-max=0.30 (overridable)
  --raw-p-threshold <float> 0.01         Raw p threshold if no FDR
  --fdr                  auto            Multiple testing: auto | bh | storey | none
  --max-prev-diff        off             Also emit “perfect” k-mers (100% in one trait group, 0% in the other)

Minimal example:
  ./BiTUGA.sh --input-dir data/reads --trait-info trait_info.tsv --out-dir run1 --threads 12 --mem-gb 16

Full option list: ./BiTUGA.sh --help+
USAGE
  exit 1
}

usage_plus() {
  cat >&2 <<'USAGE'
Usage:

  ./BiTUGA.sh --input-dir data/reads1 --trait-info trait_info.tsv [other options…]

 ./BiTUGA.sh \
    [--make-metadata] [--make-metadata-only] \
    [--input-dir <dir> ...] [--meta-output <file>] [--trait-info <file>] [--exclude-sample-id \"id1,id2\"] \
    [--metadata <file>] \
    [--out-dir <dir>] \
    [--tmpdir <dir>] \
    [--kmc-bin <dir>] \
    [--k <int>] [--ci <int>] [--cx <int>] \
    [--threads <int>] [--mem-gb <int>] \
    [--prev-min <float>] [--prev-max <float>] \
    [--group-min <float>] [--group-max <float>] \
    [--perfect-only] [--perfect-extra] \
    [--shannon-entropy-kmer <float>] \
    [--start-pass <0-5>] [--end-pass <0-5>]

Argument groups:

  Step 0 — optional metadata builder (runs only if --make-metadata is given):
    Option                Default           Description
    --make-metadata       off               Build metadata before pipeline (then continue)
    --make-metadata-only  off               Run only metadata builder and exit
    --input-dir <dir>     -                 One or more dirs to scan for FASTQ/FASTQ.GZ (repeatable)
    --meta-output <file>  ./metadata.tsv    Output metadata path
    --trait-info <file>   ./trait_info.tsv  Trait info table
    --exclude-sample-id   -                 Optional sample IDs to exclude (comma/space separated)

  Shared resources:
    Option                Default           Description
    --threads <int>       12                Max CPU threads
    --mem-gb <int>        16                Max RAM in GB

  Pipeline control:
    --start-pass <0-5>    0                 First pass to run
    --end-pass <0-5>      5                 Last pass to run

  Step I — pooled k-mers:
    --metadata <file>     ./metadata.tsv    Input metadata
    --out-dir <dir>       .                 Output and working directory
    --tmpdir <dir>        -                 Temp directory
    --k <int>             31                k-mer length
    --ci <int>            2                 KMC min counter (min enforced: 2)
    --cx <int>            4294967295        KMC max counter

  Step II — k-mer selection:
    --prev-min <float>    0.25              Global min prevalence (overrides adaptive)
    --prev-max <float>    0.75              Global max prevalence (overrides adaptive)
    --adaptive-prev auto|off auto    Adaptive defaults for Pass II (off = static 0.25/0.75)
    --group-contrast      off               Shortcut: set group-min=0.50 and group-max=0.30 (can be overridden)
    --group-min <float>   0.00              Asymmetric: one trait ≥ group-min
    --group-max <float>   1.00              Asymmetric: other trait ≤ group-max
    --max-prev-diff-only  off               Emit only “perfect” k-mers (100% vs 0%) to main FASTA
    --max-prev-diff       off               Also write all “perfect” k-mers to results/kmers_perfect_contrast.fasta
    --shannon-entropy-kmer 1.5             K-mer Shannon entropy threshold (≤0 disables)
    --homopolymer-max-frac 0.33            Max allowed homopolymer run as fraction of k (cap=15)
    --no-homopolymer-filter off            Disable homopolymer filter

  Step III — unitigs:
    --augment-singleton-unitigs off         Keep all unitigs (no length filter, min_mem=k) and append k-length unitigs from candidate_kmers

  Step IV — nt_mini_matcher:
    --nt-threads          (threads)         Threads for nt_mini_matcher (default = --threads, capped at 16)
    --min-carrier         3                 Min carriers to mark unitig present (scale with expected coverage)
    --anchor-bp           12                Anchor size for MEM stitching
    --query-split-size    262144            Max reads per batch (nt_mini_matcher default; override with care)
    --kmer-length         19                k for MEM search
    --min-mem-length      k+1               Min MEM length
    --reverse-complement  for_reference     Reverse complement mode
    --no-index-cache      off               Disable index cache
    --short-frac          1.0               Threshold for short unitigs
    --medium-frac         0.95              Threshold for medium unitigs
    --long-frac           0.90              Threshold for long unitigs
    --max-replicate-ref   -                 Limit replicated seeds on reference (pass-through to nt_mini_matcher)
    --max-replicate-qry   -                 Limit replicated seeds on query (pass-through to nt_mini_matcher)
    --short-min-bp        auto              Min bp to call short present (auto from unitig lengths; override requires setting all classes explicitly)
    --medium-min-bp       auto              Min bp to call medium present (auto from unitig lengths; override requires setting all classes explicitly)
    --long-min-bp         auto              Min bp to call long present (auto from unitig lengths; override requires setting all classes explicitly)
    --bridge-min-span     180               Min span for bridging pairs
    --short-min-segment   auto              Min segment len (short; auto from unitig lengths; override requires setting all classes explicitly)
    --medium-min-segment  auto              Min segment len (medium; auto from unitig lengths; override requires setting all classes explicitly)
    --long-min-segment    auto              Min segment len (long; auto from unitig lengths; override requires setting all classes explicitly)
    --min-bridging-pairs  2                 Min read pairs to bridge gaps
    --full-gapless-long   off               Enable gapless mode for long unitigs

  Step V — Fisher tests:
    --fisher-min-global-prev 0.10           Min global prevalence fraction to test unitig (overrides adaptive)
    --fisher-max-global-prev 0.90           Max global prevalence fraction to test unitig (overrides adaptive)
    --adaptive-fisher-prev auto|off auto    Adaptive prevalence thresholds for asymmetric groups (off = static defaults)
    --fisher-min-delta-prev  0.0            Min absolute prevalence difference to test unitig
    --fisher-min-bias-prev   0.0            Min prevalence in biased trait for FASTA output
    --fdr-alpha              0.05           Alpha for BH if used
    --fdr                    auto           Multiple testing: auto | bh | storey | none
    --raw-p-threshold        0.01           Raw p threshold when fdr-mode=none/auto fallback
    --fisher-split-by-bias   off            Write two separate FASTA (for each trait group)

  Plotting:
    --skip-plots           off            Disable all plotting
USAGE
  exit 1
}

METADATA=""
OUTDIR="${ROOT_DIR}"
USER_TMPDIR=""
K=""; CI=""; CX=""; THREADS=""; MEM_GB=""
KMC_BIN=""
PREV_MIN=""; PREV_MAX=""
ADAPTIVE_MODE="auto"
GROUP_MIN=""; GROUP_MAX=""
GROUP_CONTRAST="0"
PERFECT_ONLY="0"; PERFECT_EXTRA="0"
ENTROPY_KMER=""
HOMO_FRAC=""
HOMO_OFF="0"
ENTROPY_UNITIG=""
METADATA_SET="0"
NT_THREADS=""
P4_CARRIER=""
P4_ANCHOR=""
P4_QSPLIT=""
P4_KMERLEN=""
P4_MINMEM=""
P4_RC=""
P4_SFRAC=""
P4_MFRAC=""
P4_LFRAC=""
P4_SMINBP=""
P4_MMINBP=""
P4_LMINBP=""
P4_BRIDGE=""
P4_SMINSEG=""
P4_MMINSEG=""
P4_LMINSEG=""
P4_MINBRIDGE=""
P4_GAPLESS="0"
P4_AUG_SINGLE="0"
P4_MAXREP_REF=""
P4_MAXREP_QRY=""
FISH_MIN_GLOBAL=""
FISH_MIN_DELTA=""
FISH_MAX_GLOBAL=""
FISH_ADAPT="auto"
FISH_MIN_BIAS=""
FISH_ALPHA=""
FISH_SIG_MODE=""
FISH_PTHR="0.01"
FISH_SPLIT="0"
P4_NO_CACHE="0"
SKIP_PLOTS="0"
PLOT_OK="0"
START_PASS=0
END_PASS=5
MAKE_META_FULL="0"
MAKE_META_ONLY="0"
META_DIRS=()
META_OUT=""
META_TRAITS=""
META_EXCLUDE=()
FULL_CMD=("$0" "$@")

validate_int_min(){
  local name="$1" val="$2" min="$3"
  if ! [[ "${val}" =~ ^[0-9]+$ ]]; then
    echo "[ERR] --${name} must be an integer >=${min}" >&2
    exit 1
  fi
  if (( val < min )); then
    echo "[ERR] --${name} must be >=${min}" >&2
    exit 1
  fi
}
validate_float_01(){
  local name="$1" val="$2"
  if ! [[ "${val}" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    echo "[ERR] --${name} must be a number between 0 and 1" >&2
    exit 1
  fi
  if ! awk -v v="${val}" 'BEGIN{exit (v>=0 && v<=1)?0:1}'; then
    echo "[ERR] --${name} must be between 0 and 1" >&2
    exit 1
  fi
}
validate_float_nonneg(){
  local name="$1" val="$2"
  if ! [[ "${val}" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    echo "[ERR] --${name} must be a non-negative number" >&2
    exit 1
  fi
  if ! awk -v v="${val}" 'BEGIN{exit (v>=0)?0:1}'; then
    echo "[ERR] --${name} must be >=0" >&2
    exit 1
  fi
}
validate_int_range(){
  local name="$1" val="$2" min="$3" max="$4"
  if ! [[ "${val}" =~ ^[0-9]+$ ]]; then
    echo "[ERR] --${name} must be an integer between ${min} and ${max}" >&2
    exit 1
  fi
  if (( val < min || val > max )); then
    echo "[ERR] --${name} must be between ${min} and ${max}" >&2
    exit 1
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --help) usage ;;
    --help+) usage_plus ;;
    --make-metadata) MAKE_META_FULL="1"; shift 1;;
    --make-metadata-only) MAKE_META_ONLY="1"; shift 1;;
    --input-dir) META_DIRS+=("$2"); shift 2;;
    --meta-dir) META_DIRS+=("$2"); shift 2;;
    --meta-output) META_OUT="$2"; shift 2;;
    --trait-info) META_TRAITS="$2"; shift 2;;
    --meta-traits) META_TRAITS="$2"; shift 2;;
    --fdr-alpha) FISH_ALPHA="$2"; shift 2;;
    --raw-p-threshold) FISH_PTHR="$2"; shift 2;;
    --fdr) FISH_SIG_MODE="$2"; shift 2;;
    --sig-mode) FISH_SIG_MODE="$2"; shift 2;;
    --max-prev-diff) PERFECT_EXTRA="1"; shift 1;;
    --meta-exclude) META_EXCLUDE+=("$2"); shift 2;;
    --exclude-sample-id) META_EXCLUDE+=("$2"); shift 2;;
    --metadata) METADATA="$2"; METADATA_SET="1"; shift 2;;
    --outdir|--out-dir)   OUTDIR="$2"; shift 2;;
    --tmpdir)   USER_TMPDIR="$2"; shift 2;;
    --kmc-bin)  KMC_BIN="$2"; shift 2;;
    --k)        K="$2"; shift 2;;
    --ci)       CI="$2"; shift 2;;
    --cx)       CX="$2"; shift 2;;
    --threads)  THREADS="$2"; shift 2;;
    --mem-gb)   MEM_GB="$2"; shift 2;;
    --prev-min) PREV_MIN="$2"; shift 2;;
    --prev-max) PREV_MAX="$2"; shift 2;;
    --adaptive-prev|--adaptive-prev-thresholds) ADAPTIVE_MODE="$2"; shift 2;;
    --disable-adaptive-defaults|--force-static-defaults) ADAPTIVE_MODE="off"; shift 1;;
    --group-contrast) GROUP_CONTRAST="1"; shift 1;;
    --group-min) GROUP_MIN="$2"; shift 2;;
    --group-max) GROUP_MAX="$2"; shift 2;;
    --perfect-only)  PERFECT_ONLY="1"; shift 1;;
    --perfect-extra) PERFECT_EXTRA="1"; shift 1;;
    --max-prev-diff-only) PERFECT_ONLY="1"; shift 1;;
    --max-prev-diff) PERFECT_EXTRA="1"; shift 1;;
    --shannon-entropy-kmer) ENTROPY_KMER="$2"; shift 2;;
    --entropy-min) ENTROPY_KMER="$2"; shift 2;; # alias
    --homopolymer-max-frac) HOMO_FRAC="$2"; shift 2;;
    --no-homopolymer-filter) HOMO_OFF="1"; shift 1;;
    --max-replicate-ref) P4_MAXREP_REF="$2"; shift 2;;
    --max-replicate-qry) P4_MAXREP_QRY="$2"; shift 2;;
    --nt-threads) NT_THREADS="$2"; shift 2;;
    --carrier-threshold) P4_CARRIER="$2"; shift 2;;
    --min-carrier) P4_CARRIER="$2"; shift 2;;
    --anchor-bp) P4_ANCHOR="$2"; shift 2;;
    --query-split-size) P4_QSPLIT="$2"; shift 2;;
    --kmer-length) P4_KMERLEN="$2"; shift 2;;
    --min-mem-length) P4_MINMEM="$2"; shift 2;;
    --reverse-complement) P4_RC="$2"; shift 2;;
    --no-index-cache) P4_NO_CACHE="1"; shift 1;;
    --short-frac) P4_SFRAC="$2"; shift 2;;
    --medium-frac) P4_MFRAC="$2"; shift 2;;
    --long-frac) P4_LFRAC="$2"; shift 2;;
    --short-min-bp) P4_SMINBP="$2"; shift 2;;
    --medium-min-bp) P4_MMINBP="$2"; shift 2;;
    --long-min-bp) P4_LMINBP="$2"; shift 2;;
    --bridge-min-span) P4_BRIDGE="$2"; shift 2;;
    --short-min-segment) P4_SMINSEG="$2"; shift 2;;
    --medium-min-segment) P4_MMINSEG="$2"; shift 2;;
    --long-min-segment) P4_LMINSEG="$2"; shift 2;;
    --min-bridging-pairs) P4_MINBRIDGE="$2"; shift 2;;
    --full-gapless-long) P4_GAPLESS="1"; shift 1;;
    --augment-singleton-unitigs) P4_AUG_SINGLE="1"; shift 1;;
    --augment-k-singletons) P4_AUG_SINGLE="1"; shift 1;; # alias
    --shannon-entropy-unitig) ENTROPY_UNITIG="$2"; shift 2;;
    --fisher-min-global-prev) FISH_MIN_GLOBAL="$2"; shift 2;;
    --fisher-max-global-prev) FISH_MAX_GLOBAL="$2"; shift 2;;
    --adaptive-fisher-prev|--fisher-adaptive-prev-thresholds) FISH_ADAPT="$2"; shift 2;;
    --fisher-min-delta-prev) FISH_MIN_DELTA="$2"; shift 2;;
    --fisher-min-bias-prev) FISH_MIN_BIAS="$2"; shift 2;;
    --fisher-alpha) FISH_ALPHA="$2"; shift 2;;
    --fisher-fdr-mode) FISH_SIG_MODE="$2"; shift 2;;
    --fisher-p-threshold) FISH_PTHR="$2"; shift 2;;
    --fisher-split-by-bias) FISH_SPLIT="1"; shift 1;;
    --skip-plots) SKIP_PLOTS="1"; shift 1;;
    --start-pass) START_PASS="$2"; shift 2;;
    --end-pass) END_PASS="$2"; shift 2;;
    -h|--help)  usage;;
    *) echo "[ERR] Unknown argument: $1" >&2; usage;;
  esac
done

if [[ -z "${FISH_SIG_MODE}" ]]; then
  FISH_SIG_MODE="auto"
else
  # normalize; unknown -> error
  FISH_SIG_MODE="$(printf '%s' "${FISH_SIG_MODE}" | tr '[:upper:]' '[:lower:]')"
  case "${FISH_SIG_MODE}" in
    auto|bh|storey|none) : ;;
    *)
      echo "[ERR] --fdr must be one of: auto, bh, storey, none" >&2
      exit 1
      ;;
  esac
fi

if [[ "${GROUP_CONTRAST}" == "1" ]]; then
  [[ -z "${GROUP_MIN}" ]] && GROUP_MIN="0.50"
  [[ -z "${GROUP_MAX}" ]] && GROUP_MAX="0.30"
fi

# Argument validation
[[ -n "${THREADS}"    ]] && validate_int_min "threads" "${THREADS}" 1
[[ -n "${NT_THREADS}" ]] && validate_int_min "nt-threads" "${NT_THREADS}" 1
[[ -n "${MEM_GB}"     ]] && validate_int_min "mem-gb" "${MEM_GB}" 1
[[ -n "${P4_QSPLIT}"  ]] && validate_int_min "query-split-size" "${P4_QSPLIT}" 1
[[ -n "${K}"          ]] && validate_int_min "k" "${K}" 1
[[ -n "${CI}"         ]] && validate_int_min "ci" "${CI}" 1
[[ -n "${CX}"         ]] && validate_int_min "cx" "${CX}" 1
[[ -n "${P4_ANCHOR}"  ]] && validate_int_min "anchor-bp" "${P4_ANCHOR}" 1
[[ -n "${P4_MINBRIDGE}" ]] && validate_int_min "min-bridging-pairs" "${P4_MINBRIDGE}" 1
[[ -n "${P4_BRIDGE}"  ]] && validate_int_min "bridge-min-span" "${P4_BRIDGE}" 1
[[ -n "${P4_CARRIER}" ]] && validate_int_min "min-carrier" "${P4_CARRIER}" 1
[[ -n "${P4_KMERLEN}" ]] && validate_int_min "kmer-length" "${P4_KMERLEN}" 1
[[ -n "${P4_MINMEM}"  ]] && validate_int_min "min-mem-length" "${P4_MINMEM}" 1
[[ -n "${P4_SMINBP}"  ]] && validate_int_min "short-min-bp" "${P4_SMINBP}" 1
[[ -n "${P4_MMINBP}"  ]] && validate_int_min "medium-min-bp" "${P4_MMINBP}" 1
[[ -n "${P4_LMINBP}"  ]] && validate_int_min "long-min-bp" "${P4_LMINBP}" 1
[[ -n "${P4_SMINSEG}" ]] && validate_int_min "short-min-segment" "${P4_SMINSEG}" 1
[[ -n "${P4_MMINSEG}" ]] && validate_int_min "medium-min-segment" "${P4_MMINSEG}" 1
[[ -n "${P4_LMINSEG}" ]] && validate_int_min "long-min-segment" "${P4_LMINSEG}" 1

[[ -n "${PREV_MIN}"   ]] && validate_float_01 "prev-min" "${PREV_MIN}"
[[ -n "${PREV_MAX}"   ]] && validate_float_01 "prev-max" "${PREV_MAX}"
[[ -n "${GROUP_MIN}"  ]] && validate_float_01 "group-min" "${GROUP_MIN}"
[[ -n "${GROUP_MAX}"  ]] && validate_float_01 "group-max" "${GROUP_MAX}"
[[ -n "${FISH_MIN_GLOBAL}" ]] && validate_float_01 "fisher-min-global-prev" "${FISH_MIN_GLOBAL}"
[[ -n "${FISH_MAX_GLOBAL}" ]] && validate_float_01 "fisher-max-global-prev" "${FISH_MAX_GLOBAL}"
[[ -n "${FISH_MIN_DELTA}"  ]] && validate_float_01 "fisher-min-delta-prev" "${FISH_MIN_DELTA}"
[[ -n "${FISH_MIN_BIAS}"   ]] && validate_float_01 "fisher-min-bias-prev" "${FISH_MIN_BIAS}"
[[ -n "${FISH_PTHR}"       ]] && validate_float_01 "fisher-p-threshold" "${FISH_PTHR}"
[[ -n "${FISH_ALPHA}"      ]] && validate_float_01 "fdr-alpha" "${FISH_ALPHA}"
[[ -n "${HOMO_FRAC}"       ]] && validate_float_01 "homopolymer-max-frac" "${HOMO_FRAC}"
[[ -n "${P4_SFRAC}"        ]] && validate_float_01 "short-frac" "${P4_SFRAC}"
[[ -n "${P4_MFRAC}"        ]] && validate_float_01 "medium-frac" "${P4_MFRAC}"
[[ -n "${P4_LFRAC}"        ]] && validate_float_01 "long-frac" "${P4_LFRAC}"
[[ -n "${ENTROPY_KMER}"    ]] && validate_float_nonneg "shannon-entropy-kmer" "${ENTROPY_KMER}"
[[ -n "${ENTROPY_UNITIG}"  ]] && validate_float_nonneg "shannon-entropy-unitig" "${ENTROPY_UNITIG}"
[[ -n "${P4_KMERLEN}" && -n "${K}" ]] && { if (( P4_KMERLEN < K )); then echo "[ERR] --kmer-length must be >= --k" >&2; exit 1; fi; }
if [[ -n "${PREV_MIN}" && -n "${PREV_MAX}" ]]; then
  if awk -v a="${PREV_MIN}" -v b="${PREV_MAX}" 'BEGIN{exit (a<=b)?0:1}'; then :; else
    echo "[ERR] --prev-min must be <= --prev-max" >&2
    exit 1
  fi
fi
if [[ -n "${GROUP_MIN}" && -n "${GROUP_MAX}" ]]; then
  if awk -v a="${GROUP_MIN}" -v b="${GROUP_MAX}" 'BEGIN{exit (a<=b)?0:1}'; then :; else
    echo "[ERR] --group-min must be <= --group-max" >&2
    exit 1
  fi
fi
if [[ -n "${FISH_MIN_GLOBAL}" && -n "${FISH_MAX_GLOBAL}" ]]; then
  if awk -v a="${FISH_MIN_GLOBAL}" -v b="${FISH_MAX_GLOBAL}" 'BEGIN{exit (a<=b)?0:1}'; then :; else
    echo "[ERR] --fisher-min-global-prev must be <= --fisher-max-global-prev" >&2
    exit 1
  fi
fi
if [[ -n "${START_PASS}" ]]; then validate_int_range "start-pass" "${START_PASS}" 0 5; fi
if [[ -n "${END_PASS}" ]]; then validate_int_range "end-pass" "${END_PASS}" 0 5; fi
if [[ -n "${START_PASS}" && -n "${END_PASS}" && "${START_PASS}" -gt "${END_PASS}" ]]; then
  echo "[ERR] --start-pass cannot exceed --end-pass" >&2
  exit 1
fi

# Default metadata location: use OUTDIR/metadata.tsv if OUTDIR is set; otherwise stay in ROOT_DIR
if [[ -z "${METADATA}" ]]; then
  if [[ "${OUTDIR%/}" != "${ROOT_DIR%/}" ]]; then
    METADATA="${OUTDIR%/}/metadata.tsv"
  else
    METADATA="${ROOT_DIR}/metadata.tsv"
  fi
fi

case "${METADATA}" in /*) : ;; *) METADATA="${ROOT_DIR}/${METADATA}";; esac
case "${OUTDIR}"   in /*) : ;; *) OUTDIR="${ROOT_DIR}/${OUTDIR}";; esac

META_OUT_EFF="${META_OUT:-${METADATA}}"
META_TRAITS_EFF="${META_TRAITS:-${ROOT_DIR}/trait_info.tsv}"
case "${META_OUT_EFF}"    in /*) : ;; *) META_OUT_EFF="${ROOT_DIR}/${META_OUT_EFF}";; esac
case "${META_TRAITS_EFF}" in /*) : ;; *) META_TRAITS_EFF="${ROOT_DIR}/${META_TRAITS_EFF}";; esac
META_DIRS_ABS=()
for d in "${META_DIRS[@]}"; do
  case "${d}" in /*) META_DIRS_ABS+=("${d}") ;; *) META_DIRS_ABS+=("${ROOT_DIR}/${d}") ;; esac
done

ensure_meta(){
  [[ -x "${PASS0_BIN}" ]] || { echo "[ERR] not executable: ${PASS0_BIN}" >&2; exit 1; }
  if [[ ${#META_DIRS_ABS[@]} -eq 0 ]]; then
    echo "[ERR] --meta-dir <dir> required to build metadata (Pass 0)" >&2
    exit 1
  fi
  if [[ ! -f "${META_TRAITS_EFF}" && -z "${META_TRAITS:-}" ]]; then
    for cand in "${ROOT_DIR}/trait_info.tsv" "${ROOT_DIR}/trait_info.txt" "${ROOT_DIR}/trait_info.csv" "${ROOT_DIR}/trait_info"; do
      if [[ -f "${cand}" ]]; then META_TRAITS_EFF="${cand}"; break; fi
    done
  fi
  if [[ ! -f "${META_TRAITS_EFF}" ]]; then
    echo "[ERR] trait info not found (searched: ${META_TRAITS_EFF})" >&2
    exit 1
  fi
  # sanity: at least two distinct trait values
  if ! python3 - "$META_TRAITS_EFF" <<'PY'
import sys, re
path=sys.argv[1]
splitter=re.compile(r"[,;\t]+|\s+")
traits=set()
with open(path,"r",encoding="utf-8") as fh:
    header_seen=False
    for line in fh:
        line=line.strip()
        if not line or line.startswith("#"): continue
        parts=[p for p in splitter.split(line) if p]
        if not parts: continue
        if not header_seen and parts[0].lower() in {"id","sample_id","sample"}:
            header_seen=True
            continue
        if len(parts)>=2:
            traits.add(parts[1].strip().lower())
if len(traits) < 2:
    sys.stderr.write(f"[ERR] trait info must contain at least two distinct trait values: {path}\n")
    sys.exit(2)
PY
  then
    echo "[ERR] trait_info sanity check failed" >&2
    exit $?
  fi
  CMD0=( python3 "${PASS0_BIN}" "-o" "${META_OUT_EFF}" "-t" "${META_TRAITS_EFF}" )
  for d in "${META_DIRS_ABS[@]}"; do CMD0+=( "-d" "${d}" ); done
  for ex in "${META_EXCLUDE[@]:-}"; do CMD0+=( "--exclude" "${ex}" ); done
  set +e
  { "${CMD0[@]}" 2>&1 | while IFS= read -r line; do logi "$line"; done; }
  rc0=${PIPESTATUS[0]}
  set -e
  if [[ ${rc0} -ne 0 ]]; then
    echo "[ERR] Pass 0 (make metadata) failed rc=${rc0}" >&2
    exit "${rc0}"
  fi
  METADATA="${META_OUT_EFF}"
}

LOGROOT="${OUTDIR%/}/logging"
MAINLOG="${LOGROOT}/main.log"
META_CACHE="${LOGROOT}/.metadata_path"
mkdir -p "${LOGROOT}"

logi_cmd(){
  local msg="$*"
  local line="[$(ts)] [INFO] ${msg}"
  echo "${line}" | tee -a "${MAINLOG}"
}

ts(){ date '+%Y-%m-%d %H:%M:%S'; }
rel_out(){
  local msg="$*"
  local base="${OUTDIR%/}"
  if [[ -n "${base}" ]]; then
    msg="${msg//${base}\//./}"
    msg="${msg//${base}/.}"
  fi
  local root="${ROOT_DIR%/}"
  if [[ -n "${root}" ]]; then
    msg="${msg//${root}\//./}"
    msg="${msg//${root}/.}"
  fi
  printf '%s' "${msg}"
}
logi(){
  local msg
  msg="$(rel_out "$*")"
  line="[$(ts)] [INFO] ${msg}"
  echo "${line}" | tee -a "${MAINLOG}"
}
logw(){
  local msg
  msg="$(rel_out "$*")"
  line="[$(ts)] [WARN] ${msg}"
  echo "${line}" | tee -a "${MAINLOG}"
}
loge(){
  local msg
  msg="$(rel_out "$*")"
  line="[$(ts)] [ERROR] ${msg}"
  echo "${line}" | tee -a "${MAINLOG}" >&2
}
logi_cmd "Command: ${FULL_CMD[*]}"
mark_path(){ echo "${LOGROOT:-.}/.pass$1.done"; }
clean_pass(){
  local p="$1"
  case "$p" in
    1)
      logi "Cleaning outputs for Pass I"
      rm -rf "${OUTDIR%/}/results/pooled_kmc_files" "${OUTDIR%/}/results/statistics/kmer_stats_per_sample.tsv" \
             "${OUTDIR%/}/results/statistics/kmer_totals.tsv" "${LOGROOT}/build_pooled_kmc_files.log" 2>/dev/null || true
      ;;
    2)
      logi "Cleaning outputs for Pass II"
      rm -f "${OUTDIR%/}/results/candidate_kmers.fasta" "${OUTDIR%/}/results/candidate_kmers.fasta.gz" \
            "${OUTDIR%/}/results/statistics/group_prevalence_histogram.txt" \
            "${LOGROOT}/candidate_kmers.log" 2>/dev/null || true
      ;;
    3)
      logi "Cleaning outputs for Pass III"
      rm -f "${OUTDIR%/}/results/unitigs.fa" "${OUTDIR%/}/results/unitigs.fasta" \
            "${OUTDIR%/}/results/statistics/unitig_statistics.txt" \
            "${LOGROOT}/build_unitigs.log" 2>/dev/null || true
      rm -rf "${OUTDIR%/}/results/unitigs_build_"* 2>/dev/null || true
      ;;
    4)
      logi "Cleaning outputs for Pass IV"
      rm -rf "${OUTDIR%/}/results/unitig_matcher" 2>/dev/null || true
      rm -f "${OUTDIR%/}/results/unitig_presence.tsv" "${OUTDIR%/}/results/statistics/unitig_matches.txt" \
            "${LOGROOT}/unitig_matching.log" 2>/dev/null || true
      ;;
    5)
      logi "Cleaning outputs for Pass V"
      rm -f "${OUTDIR%/}/results/statistics/unitig_association_stats.txt" \
            "${OUTDIR%/}/results/statistics/unitig_association_details.tsv" \
            "${OUTDIR%/}/results/statistics/unitig_fishers.txt" \
            "${OUTDIR%/}/results/statistics/unitigs_fishers_detailed.tsv" \
            "${OUTDIR%/}/results/significant_unitigs.fa" \
            "${OUTDIR%/}/results/significant_unitigs_trait"*".fa" 2>/dev/null || true
      rm -rf "${OUTDIR%/}/results/statistics/plots" 2>/dev/null || true
      ;;
  esac
  rm -f "$(mark_path "$p")" 2>/dev/null || true
}
outputs_exist(){
  local p="$1"
  case "$p" in
    1) [[ -d "${OUTDIR%/}/results/pooled_kmc_files" ]] && return 0 ;;
    2) [[ -f "${OUTDIR%/}/results/candidate_kmers.fasta" || -f "${OUTDIR%/}/results/candidate_kmers.fasta.gz" ]] && return 0 ;;
    3) [[ -f "${OUTDIR%/}/results/unitigs.fa" || -f "${OUTDIR%/}/results/unitigs.fasta" ]] && return 0 ;;
    4) [[ -f "${OUTDIR%/}/results/unitig_presence.tsv" ]] && return 0 ;;
    5) [[ -f "${OUTDIR%/}/results/statistics/unitig_association_stats.txt" ]] && return 0 ;;
  esac
  return 1
}

validate_range(){
  local v="$1"; [[ "$v" =~ ^[0-9]+$ ]] || { loge "start-pass/end-pass must be integers"; exit 1; }
  [[ "$v" -ge 0 && "$v" -le 5 ]] || { loge "start-pass/end-pass must be between 0 and 5"; exit 1; }
}
validate_range "${START_PASS}"
validate_range "${END_PASS}"
if [[ "${START_PASS}" -gt "${END_PASS}" ]]; then
  loge "start-pass (${START_PASS}) cannot be greater than end-pass (${END_PASS})"
  exit 1
fi

should_run_pass(){ local p="$1"; [[ "$p" -ge "${START_PASS}" && "$p" -le "${END_PASS}" ]]; }

if [[ "${START_PASS}" -gt 0 ]]; then
  if [[ "${METADATA_SET}" != "1" && ! -f "${METADATA}" && -f "${META_CACHE}" ]]; then
    cached_meta="$(cat "${META_CACHE}" 2>/dev/null || true)"
    if [[ -n "${cached_meta}" && -f "${cached_meta}" ]]; then
      logi "start-pass>=1: metadata not found at ${METADATA}, using cached path ${cached_meta}"
      METADATA="${cached_meta}"
    fi
  fi
  [[ -f "${METADATA}" ]] || { loge "start-pass>=1 requires metadata at ${METADATA}"; exit 1; }
  for ((p=1; p<START_PASS; p++)); do
    mk="$(mark_path "$p")"
    if [[ ! -f "${mk}" ]]; then
      loge "start-pass=${START_PASS} requested but checkpoint missing for pass ${p}; run earlier passes first"
      exit 1
    fi
  done
  for ((p=START_PASS; p<=5; p++)); do
    pass_mark="$(mark_path "$p")"
    if [[ -f "${pass_mark}" ]]; then
      logi "start-pass=${START_PASS}: cleaning checkpoint for pass ${p}"
      clean_pass "$p"
      continue
    fi
    if outputs_exist "$p"; then
      logi "start-pass=${START_PASS}: cleaning outputs for pass ${p}"
      clean_pass "$p"
    fi
  done
fi

if should_run_pass 0; then
  if [[ "${MAKE_META_ONLY}" == "1" ]]; then
    if [[ -f "${METADATA}" && "${MAKE_META_FULL}" != "1" ]]; then
      logi "Pass 0 skipped (metadata exists)"
      touch "$(mark_path 0)"
    else
      logi "Pass 0 === MAKE METADATA ==="
      ensure_meta
      logi "Pass 0 end rc=0 (metadata=$(basename "${METADATA}"))"
      mkdir -p "$(dirname "${META_CACHE}")"
      printf "%s\n" "${METADATA}" > "${META_CACHE}"
      touch "$(mark_path 0)"
    fi
    exit 0
  fi

  if [[ ! -f "${METADATA}" || "${MAKE_META_FULL}" == "1" ]]; then
    logi "Pass 0 === MAKE METADATA ==="
    ensure_meta
    logi "Pass 0 end rc=0 (metadata=$(basename "${METADATA}"))"
    mkdir -p "$(dirname "${META_CACHE}")"
    printf "%s\n" "${METADATA}" > "${META_CACHE}"
    touch "$(mark_path 0)"
  fi
fi

if [[ -f "${METADATA}" ]]; then
  mkdir -p "$(dirname "${META_CACHE}")"
  printf "%s\n" "${METADATA}" > "${META_CACHE}"
fi

[[ -f "${METADATA}" ]]     || { loge "metadata not found: ${METADATA}"; exit 1; }

[[ -x "${BUILD}" ]]        || { echo "[ERR] not executable: ${BUILD}" >&2; exit 1; }
[[ -f "${METADATA}" ]]     || { echo "[ERR] metadata not found: ${METADATA}" >&2; exit 1; }
[[ -x "${PASS2_BIN}" ]]    || { echo "[ERR] not executable: ${PASS2_BIN}" >&2; exit 1; }
[[ -x "${NM}" ]]           || { echo "[ERR] not executable: ${NM}" >&2; exit 1; }
[[ -x "${PERF}" ]]         || { echo "[ERR] not executable: ${PERF}" >&2; exit 1; }
[[ -x "${PASS3_PY}" ]]     || { echo "[ERR] not executable: ${PASS3_PY}" >&2; exit 1; }
[[ -x "${PASS4_BIN}" ]]    || { echo "[ERR] not executable: ${PASS4_BIN}" >&2; exit 1; }
[[ -x "${PASS5_BIN}" ]]    || { echo "[ERR] not executable: ${PASS5_BIN}" >&2; exit 1; }

ul_soft="$(ulimit -Sn || echo 256)"; if [[ "${ul_soft}" -lt 4096 ]]; then ulimit -n 4096 2>/dev/null || true; fi

CLEAN_LINK="false"
TMP_TARGET=""
if [[ -n "${USER_TMPDIR}" ]]; then
  case "${USER_TMPDIR}" in /*) : ;; *) USER_TMPDIR="${ROOT_DIR}/${USER_TMPDIR}";; esac
  mkdir -p "${USER_TMPDIR}"
  ts_tag="$(date '+%Y%m%d_%H%M%S')"
  TMP_TARGET="${USER_TMPDIR%/}/sexfinder_tmp_${ts_tag}_$$"
  mkdir -p "${TMP_TARGET}"
  if [[ -L "${OUTDIR}/tmp" ]]; then
    rm -f "${OUTDIR}/tmp"
  elif [[ -d "${OUTDIR}/tmp" ]]; then
    [[ -z "$(ls -A "${OUTDIR}/tmp")" ]] && rmdir "${OUTDIR}/tmp" || { loge "${OUTDIR}/tmp not empty"; exit 2; }
  elif [[ -e "${OUTDIR}/tmp" ]]; then
    loge "${OUTDIR}/tmp exists and is not dir/symlink"; exit 2
  fi
  ln -s "${TMP_TARGET}" "${OUTDIR}/tmp"; CLEAN_LINK="true"
fi

K_EFF="${K:-31}"
CI_EFF="${CI:-2}"
CX_EFF="${CX:-4294967295}"
if ! [[ "${CI_EFF}" =~ ^[0-9]+$ ]]; then
  loge "--ci must be an integer"; exit 1;
fi
if [[ "${CI_EFF}" -lt 2 ]]; then
  logi "--ci requested <2; clamping to 2"
  CI_EFF=2
fi
CS_EFF="2"
THREADS_EFF="${THREADS:-12}"
MEM_EFF="${MEM_GB:-16}"
if [[ "${P4_AUG_SINGLE}" == "1" ]]; then
  UNITIG_MINLEN_EFF="${K_EFF}"
else
  UNITIG_MINLEN_EFF="$((K_EFF + 1))"
fi
ENTROPY_UNITIG_EFF="${ENTROPY_UNITIG:-1.3}"
if [[ -n "${NT_THREADS}" ]]; then
  NT_THREADS_EFF="${NT_THREADS}"
else
  NT_THREADS_EFF="${THREADS_EFF}"
  if [[ "${NT_THREADS_EFF}" -gt 12 ]]; then NT_THREADS_EFF=12; fi
fi
P4_QSPLIT_EFF="${P4_QSPLIT:-${QSPLIT_DEFAULT}}"
P4_MINMEM_EFF="${P4_MINMEM:-$((K_EFF + 1))}"
if [[ -z "${P4_MINMEM:-}" && "${P4_AUG_SINGLE}" == "1" ]]; then
  P4_MINMEM_EFF="${K_EFF}"
fi

ARGS_A=( "--metadata" "${METADATA}" "--outdir" "${OUTDIR}" )
[[ -n "${K}"       ]] && ARGS_A+=( "--k" "${K}" )
ARGS_A+=( "--ci" "${CI_EFF}" )
[[ -n "${CX}"      ]] && ARGS_A+=( "--cx" "${CX}" )
ARGS_A+=( "--cs" "${CS_EFF}" )
[[ -n "${THREADS}" ]] && ARGS_A+=( "--threads" "${THREADS}" )
[[ -n "${MEM_GB}"  ]] && ARGS_A+=( "--mem-gb" "${MEM_GB}" )
[[ -n "${KMC_BIN}" ]] && ARGS_A+=( "--kmc-bin" "${KMC_BIN}" )

PASS1_MARK="$(mark_path 1)"
if ! should_run_pass 1; then
  logi "Pass I skipped (outside start/end window)"
elif [[ -f "${PASS1_MARK}" ]]; then
  logi "Pass I skipped (checkpoint exists: ${PASS1_MARK})"
else
  if outputs_exist 1; then
    logi "Pass I outputs found without checkpoint — cleaning and rerunning"
    clean_pass 1
  fi
      logi "Pass I === KMC K-MER AGGREGATION ==="
  set +e
  "${PERF}" --label "Pass I (build pooled kmc files)" -- "${BUILD}" "${ARGS_A[@]}" 2>&1 | tee -a "${MAINLOG}"
  rc_a=${PIPESTATUS[0]}
  set -e
  logi "Pass I end rc=${rc_a}"
  if [[ ${rc_a} -ne 0 ]]; then
    [[ "${CLEAN_LINK}" == "true" ]] && [[ -L "${OUTDIR}/tmp" ]] && rm -f "${OUTDIR}/tmp"
    [[ -n "${TMP_TARGET}" && -d "${TMP_TARGET}" ]] && rm -rf "${TMP_TARGET}" || true
    exit "${rc_a}"
  fi
  touch "${PASS1_MARK}"
fi

readout="$(python3 "${NM}" --in "${METADATA}" --summary)" || { loge "Failed to summarize traits from ${METADATA}"; exit 9; }
T0="$(echo "${readout}" | awk -F'\t' 'NR==1{print $1}')"
N0="$(echo "${readout}" | awk -F'\t' 'NR==1{print $2}')"
T1="$(echo "${readout}" | awk -F'\t' 'NR==1{print $3}')"
N1="$(echo "${readout}" | awk -F'\t' 'NR==1{print $4}')"

RES_DIR="${OUTDIR%/}/results"
LOG_DIR="${OUTDIR%/}/logging"
PKF_DIR="${RES_DIR}/pooled_kmc_files"
INDEX0="${PKF_DIR}/pooled_${T0}"
INDEX1="${PKF_DIR}/pooled_${T1}"
OUT_PREFIX="${RES_DIR}/candidate_kmers"
CAND_LOG="${LOG_DIR}/candidate_kmers.log"
mkdir -p "${RES_DIR}" "${LOG_DIR}"

CMD_B=( "${PASS2_BIN}"
  "--n0" "${N0}" "--n1" "${N1}"
  "--trait0" "${T0}" "--trait1" "${T1}"
  "--output" "${OUTDIR}"
  "--results-dir" "${RES_DIR}"
  "--logging-dir" "${LOG_DIR}"
  "--index0" "${INDEX0}" "--index1" "${INDEX1}"
  "--out-prefix" "${OUT_PREFIX}"
)
[[ -n "${PREV_MIN}"   ]] && CMD_B+=( "--prev-min" "${PREV_MIN}" )
[[ -n "${PREV_MAX}"   ]] && CMD_B+=( "--prev-max" "${PREV_MAX}" )
[[ -n "${ADAPTIVE_MODE}" ]] && CMD_B+=( "--adaptive-prev-thresholds" "${ADAPTIVE_MODE}" )
[[ -n "${GROUP_MIN}"  ]] && CMD_B+=( "--group-min" "${GROUP_MIN}" )
[[ -n "${GROUP_MAX}"  ]] && CMD_B+=( "--group-max" "${GROUP_MAX}" )
[[ "${PERFECT_ONLY}"  == "1" ]] && CMD_B+=( "--perfect-only" )
[[ "${PERFECT_EXTRA}" == "1" ]] && CMD_B+=( "--perfect-extra" )
[[ -n "${ENTROPY_KMER}" ]] && CMD_B+=( "--shannon-entropy-kmer" "${ENTROPY_KMER}" )
[[ "${HOMO_OFF}" == "1" ]] && CMD_B+=( "--no-homopolymer-filter" )
[[ -n "${HOMO_FRAC}" ]] && CMD_B+=( "--homopolymer-max-frac" "${HOMO_FRAC}" )

printf -v CMD_STR '%q ' "${CMD_B[@]}"

PASS2_MARK="$(mark_path 2)"
if ! should_run_pass 2; then
  logi "Pass II skipped (outside start/end window)"
elif [[ -f "${PASS2_MARK}" ]]; then
  logi "Pass II skipped (checkpoint exists: ${PASS2_MARK})"
else
  if [[ ! -f "${PASS2_MARK}" ]] && outputs_exist 2; then
    logi "Pass II outputs found without checkpoint — cleaning and rerunning"
    clean_pass 2
  fi
  logi "Pass II === K-MER FILTERING ==="
  CMD_STR_LOG="$(rel_out "${CMD_STR}")"
  logi "Pass II run: ${CMD_STR_LOG}"
  set +e
  PERF_CMD="set -o pipefail; echo \"[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] Pass II run: ${CMD_STR_LOG}\" | tee -a '${CAND_LOG}'; ${CMD_STR} 2>&1 | tee -a '${CAND_LOG}'; rc_main=\${PIPESTATUS[0]}; if [[ \${rc_main} -eq 0 && -f '${CAND_LOG}' ]]; then mkdir -p '${RES_DIR}/statistics' && cp -p '${CAND_LOG}' '${RES_DIR}/statistics/kmer_candidates.txt'; fi; exit \${rc_main}"
  "${PERF}" --label "Pass II (find candidate kmers)" -- bash -lc "${PERF_CMD}" 2>&1 | tee -a "${MAINLOG}"
  rc_b=${PIPESTATUS[0]}
  set -e
  logi "Pass II end rc=${rc_b}"
  if [[ ${rc_b} -ne 0 ]]; then
    [[ "${CLEAN_LINK}" == "true" ]] && [[ -L "${OUTDIR}/tmp" ]] && rm -f "${OUTDIR}/tmp"
    [[ -n "${TMP_TARGET}" && -d "${TMP_TARGET}" ]] && rm -rf "${TMP_TARGET}" || true
    exit "${rc_b}"
  fi
  touch "${PASS2_MARK}"
fi

PASS3_MARK="$(mark_path 3)"
if ! should_run_pass 3; then
  logi "Pass III skipped (outside start/end window)"
elif [[ -f "${PASS3_MARK}" ]]; then
  logi "Pass III skipped (checkpoint exists: ${PASS3_MARK})"
else
  if [[ ! -f "${PASS3_MARK}" ]] && outputs_exist 3; then
    logi "Pass III outputs found without checkpoint — cleaning and rerunning"
    clean_pass 3
  fi
  logi "Pass III === BCALM UNITIG ASSEMBLY ==="
  CMD_C=( python3 "${PASS3_PY}" --base-dir "${OUTDIR}" --k "${K_EFF}" --threads "${THREADS_EFF}" --max-memory-gb "${MEM_EFF}" )
  printf -v CMD3_STR '%q ' "${CMD_C[@]}"
  CMD_COMBINED="${CMD3_STR}"
  if [[ "${P4_AUG_SINGLE}" == "1" ]]; then
    logi "Pass III filter skipped (augment-singleton-unitigs enabled -> keep all unitigs)"
  else
    : # Filter deaktiviert: einfach Unitigs bauen
  fi
  logi "Pass III run: ${CMD_COMBINED}"
  set +e
  "${PERF}" --label "Pass III (build unitigs)" -- bash -lc "${CMD_COMBINED}" 2>&1 | tee -a "${MAINLOG}"
  rc_c=${PIPESTATUS[0]}
  set -e
  logi "Pass III end rc=${rc_c}"
  if [[ ${rc_c} -ne 0 ]]; then
    [[ "${CLEAN_LINK}" == "true" ]] && [[ -L "${OUTDIR}/tmp" ]] && rm -f "${OUTDIR}/tmp"
    [[ -n "${TMP_TARGET}" && -d "${TMP_TARGET}" ]] && rm -rf "${TMP_TARGET}" || true
    exit "${rc_c}"
  fi
  for tmpdir in "${OUTDIR%/}/results"/unitigs_build_*; do
    if [[ -d "${tmpdir}" ]]; then
      logi "Pass III cleanup: removing ${tmpdir}"
      rm -rf "${tmpdir}" || true
    fi
  done
  touch "${PASS3_MARK}"
fi

CMD_P4=( python3 "${PASS4_BIN}" --outdir "${OUTDIR}" --metadata "${METADATA}" --threads "${NT_THREADS_EFF}" )
[[ -n "${P4_CARRIER}"    ]] && CMD_P4+=( "--carrier-threshold" "${P4_CARRIER}" )
[[ -n "${P4_ANCHOR}"     ]] && CMD_P4+=( "--anchor-bp" "${P4_ANCHOR}" )
CMD_P4+=( "--query-split-size" "${P4_QSPLIT_EFF}" )
[[ -n "${P4_KMERLEN}"    ]] && CMD_P4+=( "--kmer-length" "${P4_KMERLEN}" )
if [[ -n "${P4_MINMEM}" ]]; then
  CMD_P4+=( "--min-mem-length" "${P4_MINMEM}" )
else
  CMD_P4+=( "--min-mem-length" "${P4_MINMEM_EFF}" )
fi
[[ -n "${P4_RC}"         ]] && CMD_P4+=( "--reverse-complement" "${P4_RC}" )
[[ -n "${P4_SFRAC}"      ]] && CMD_P4+=( "--short-frac" "${P4_SFRAC}" )
[[ -n "${P4_MFRAC}"      ]] && CMD_P4+=( "--medium-frac" "${P4_MFRAC}" )
[[ -n "${P4_LFRAC}"      ]] && CMD_P4+=( "--long-frac" "${P4_LFRAC}" )
[[ -n "${P4_SMINBP}"     ]] && CMD_P4+=( "--short-min-bp" "${P4_SMINBP}" )
[[ -n "${P4_MMINBP}"     ]] && CMD_P4+=( "--medium-min-bp" "${P4_MMINBP}" )
[[ -n "${P4_LMINBP}"     ]] && CMD_P4+=( "--long-min-bp" "${P4_LMINBP}" )
[[ -n "${P4_BRIDGE}"     ]] && CMD_P4+=( "--bridge-min-span" "${P4_BRIDGE}" )
[[ -n "${P4_SMINSEG}"    ]] && CMD_P4+=( "--short-min-segment" "${P4_SMINSEG}" )
[[ -n "${P4_MMINSEG}"    ]] && CMD_P4+=( "--medium-min-segment" "${P4_MMINSEG}" )
[[ -n "${P4_LMINSEG}"    ]] && CMD_P4+=( "--long-min-segment" "${P4_LMINSEG}" )
[[ -n "${P4_MINBRIDGE}"  ]] && CMD_P4+=( "--min-bridging-pairs" "${P4_MINBRIDGE}" )
[[ "${P4_GAPLESS}" == "1" ]] && CMD_P4+=( "--full-gapless-long" )
[[ "${P4_NO_CACHE}" == "1" ]] && CMD_P4+=( "--no-index-cache" )
[[ "${P4_AUG_SINGLE}" == "1" ]] && CMD_P4+=( "--augment-k-singletons" )
[[ -n "${P4_MAXREP_REF}" ]] && CMD_P4+=( "--max-replicate-ref" "${P4_MAXREP_REF}" )
[[ -n "${P4_MAXREP_QRY}" ]] && CMD_P4+=( "--max-replicate-qry" "${P4_MAXREP_QRY}" )

printf -v CMD_P4_STR '%q ' "${CMD_P4[@]}"
PASS4_MARK="$(mark_path 4)"
if ! should_run_pass 4; then
  logi "Pass IV skipped (outside start/end window)"
elif [[ -f "${PASS4_MARK}" ]]; then
  logi "Pass IV skipped (checkpoint exists: ${PASS4_MARK})"
else
  if [[ ! -f "${PASS4_MARK}" ]] && outputs_exist 4; then
    logi "Pass IV outputs found without checkpoint — cleaning and rerunning"
    clean_pass 4
  fi
  logi "Pass IV === UNITIG VERIFICATION ==="
  logi "Pass IV run: ${CMD_P4[0]} ${CMD_P4[1]##*/} …"
  set +e
  "${PERF}" --label "Pass IV (unitig read matching)" -- bash -lc "${CMD_P4_STR}" 2>&1 | tee -a "${MAINLOG}"
  rc_p4=${PIPESTATUS[0]}
  set -e
  logi "Pass IV end rc=${rc_p4}"
  if [[ ${rc_p4} -ne 0 ]]; then
    loge "Pass IV failed (rc=${rc_p4}). Possible OOM/kill. For large FASTQs try a smaller --query-split-size (e.g., 65536 or 131072) and cap seeds via --max-replicate-ref 100 --max-replicate-qry 100."
    loge "Retry from Pass IV with: [your BiTUGA.sh call] --query-split-size 65536 --max-replicate-ref 100 --max-replicate-qry 100 --start-pass 4"
    [[ "${CLEAN_LINK}" == "true" ]] && [[ -L "${OUTDIR}/tmp" ]] && rm -f "${OUTDIR}/tmp"
    [[ -n "${TMP_TARGET}" && -d "${TMP_TARGET}" ]] && rm -rf "${TMP_TARGET}" || true
    exit "${rc_p4}"
  fi
  touch "${PASS4_MARK}"
fi

if [[ "${CLEAN_LINK}" == "true" ]]; then
  [[ -L "${OUTDIR}/tmp" ]] && rm -f "${OUTDIR}/tmp"
  [[ -n "${TMP_TARGET}" && -d "${TMP_TARGET}" ]] && rm -rf "${TMP_TARGET}" || true
fi

CMD_P5=( python3 "${PASS5_BIN}" --outdir "${OUTDIR}" --metadata "${METADATA}" )
[[ -n "${FISH_MIN_GLOBAL}" ]] && CMD_P5+=( "--min-global-prev" "${FISH_MIN_GLOBAL}" )
[[ -n "${FISH_MAX_GLOBAL}" ]] && CMD_P5+=( "--max-global-prev" "${FISH_MAX_GLOBAL}" )
[[ -n "${FISH_ADAPT}"      ]] && CMD_P5+=( "--adaptive-prev-thresholds" "${FISH_ADAPT}" )
[[ -n "${FISH_MIN_DELTA}"  ]] && CMD_P5+=( "--min-delta-prev" "${FISH_MIN_DELTA}" )
[[ -n "${FISH_MIN_BIAS}"   ]] && CMD_P5+=( "--min-bias-prev" "${FISH_MIN_BIAS}" )
[[ -n "${FISH_ALPHA}"      ]] && CMD_P5+=( "--alpha" "${FISH_ALPHA}" )
[[ -n "${FISH_SIG_MODE}"   ]] && CMD_P5+=( "--fdr" "${FISH_SIG_MODE}" )
[[ -n "${FISH_PTHR}"       ]] && CMD_P5+=( "--p-threshold" "${FISH_PTHR}" )
[[ "${FISH_SPLIT}" == "1" ]] && CMD_P5+=( "--split-by-bias" )

printf -v CMD5_STR '%q ' "${CMD_P5[@]}"
PASS5_MARK="$(mark_path 5)"
if ! should_run_pass 5; then
  logi "Pass V skipped (outside start/end window)"
elif [[ -f "${PASS5_MARK}" ]]; then
  logi "Pass V skipped (checkpoint exists: ${PASS5_MARK})"
else
  if [[ ! -f "${PASS5_MARK}" ]] && outputs_exist 5; then
    logi "Pass V outputs found without checkpoint — cleaning and rerunning"
    clean_pass 5
  fi
  logi "Pass V start"
  logi "Pass V run: ${CMD_P5[0]} ${CMD_P5[1]##*/} …"
  set +e
  "${PERF}" --label "Pass V (fisher tests)" -- bash -lc "${CMD5_STR}" 2>&1 | tee -a "${MAINLOG}"
  rc_p5=${PIPESTATUS[0]}
  set -e
  logi "Pass V end rc=${rc_p5}"
  P5_RC="${rc_p5}"
  touch "${PASS5_MARK}"
  PLOT_OK="1"
fi

if [[ "${SKIP_PLOTS}" == "1" ]]; then
  logi "Plotting skipped (--skip-plots)"
elif [[ "${PLOT_OK}" != "1" ]]; then
  logi "Plotting skipped (Pass V not executed in this run)"
else
  if ! python3 - <<'PY'
try:
    import matplotlib  # noqa: F401
except ImportError:
    raise SystemExit(1)
PY
  then
    logw "matplotlib not installed; skipping all plotting steps (install via pip install matplotlib)"
  else
  PLOT_CMDS=()
  PLOT_MPLCONFIG="${OUTDIR%/}/results/statistics/plots/mplconfig"
  FUNNEL_SCRIPT="${ROOT_DIR}/bin/plot_pipeline_funnel.py"
  FUNNEL_PREFIX="${OUTDIR%/}/results/statistics/plots/funnel/pipeline_funnel"
  FUNNEL_CMD=( "env" "MPLCONFIGDIR=${PLOT_MPLCONFIG}" "python3" "${FUNNEL_SCRIPT}" "--run-root" "${OUTDIR}" "--out-prefix" "${FUNNEL_PREFIX}" )
  if [[ -n "${FISH_ALPHA}" ]]; then FUNNEL_CMD+=( "--alpha" "${FISH_ALPHA}" ); fi
  if [[ -n "${FISH_PTHR}" ]]; then FUNNEL_CMD+=( "--p-threshold" "${FISH_PTHR}" ); fi
  PLOT_CMDS+=( "$(printf '%q ' "${FUNNEL_CMD[@]}")" )

  KSTAT_SCRIPT="${ROOT_DIR}/bin/plot_kmer_stats.py"
  KSTAT_PREFIX="${OUTDIR%/}/results/statistics/plots/kmer_stats/kmer_stats"
  KSTAT_CMD=( "env" "MPLCONFIGDIR=${PLOT_MPLCONFIG}" "python3" "${KSTAT_SCRIPT}" "--run-root" "${OUTDIR}" "--out-prefix" "${KSTAT_PREFIX}" )
  PLOT_CMDS+=( "$(printf '%q ' "${KSTAT_CMD[@]}")" )

  KPREV_SCRIPT="${ROOT_DIR}/bin/plot_kmer_prevalence.py"
  KPREV_PREFIX="${OUTDIR%/}/results/statistics/plots/kmer_stats/kmer_group_prevalence"
  KPREV_CMD=( "env" "MPLCONFIGDIR=${PLOT_MPLCONFIG}" "python3" "${KPREV_SCRIPT}" "--run-root" "${OUTDIR}" "--out-prefix" "${KPREV_PREFIX}" )
  PLOT_CMDS+=( "$(printf '%q ' "${KPREV_CMD[@]}")" )

  PVAL_SCRIPT="${ROOT_DIR}/bin/plot_pvalue_hist.py"
  PVAL_PREFIX="${OUTDIR%/}/results/statistics/plots/pvalues/pvalues"
  PVAL_CMD=( "env" "MPLCONFIGDIR=${PLOT_MPLCONFIG}" "python3" "${PVAL_SCRIPT}" "--run-root" "${OUTDIR}" "--out-prefix" "${PVAL_PREFIX}" )
  PLOT_CMDS+=( "$(printf '%q ' "${PVAL_CMD[@]}")" )

  ULEN_SCRIPT="${ROOT_DIR}/bin/plot_unitig_lengths.py"
  ULEN_PREFIX="${OUTDIR%/}/results/statistics/plots/unitigs/unitig_lengths"
  ULEN_CMD=( "env" "MPLCONFIGDIR=${PLOT_MPLCONFIG}" "python3" "${ULEN_SCRIPT}" "--run-root" "${OUTDIR}" "--out-prefix" "${ULEN_PREFIX}" )
  PLOT_CMDS+=( "$(printf '%q ' "${ULEN_CMD[@]}")" )

  logi "Plotting start"
  set +e
  plot_shell=$(IFS='; '; echo "${PLOT_CMDS[*]}")
  "${PERF}" --label "Plotting" -- bash -lc "${plot_shell}" >> "${MAINLOG}" 2>&1
  rc_plot=${PIPESTATUS[0]}
  set -e
  if [[ ${rc_plot} -ne 0 ]]; then
    loge "Plotting failed rc=${rc_plot}"
  else
    logi "Plotting done"
    FUNNEL_TSV="${FUNNEL_PREFIX}.tsv"
    if [[ -f "${FUNNEL_TSV}" ]]; then
      while IFS= read -r line; do
        [[ -z "${line}" ]] && continue
        if [[ "${line}" == stage*$'\t'* ]]; then continue; fi
        stage_lbl=$(echo "${line}" | awk -F'\t' '{print $1}')
        stage_cnt=$(echo "${line}" | awk -F'\t' '{print $2}')
        logi "Funnel stage: ${stage_lbl} => ${stage_cnt}"
      done < "${FUNNEL_TSV}"
    fi
  fi
  if [[ -n "${PLOT_MPLCONFIG:-}" && -d "${PLOT_MPLCONFIG}" ]]; then
    rm -rf "${PLOT_MPLCONFIG}"
  fi
  fi
fi

if [[ -n "${P5_RC:-}" && "${P5_RC}" -ne 0 ]]; then
  exit "${P5_RC}"
fi

if [[ -f "${MAINLOG}" ]]; then
  logi "Perf summary aggregation start"
  if ! python3 "${ROOT_DIR}/bin/aggregate_perf_summary.py" --log "${MAINLOG}" --out "${LOGROOT}/perf_summary.txt" | tee -a "${MAINLOG}"; then
    loge "Perf summary aggregation failed"
  else
    logi "Perf summary aggregation done (written to ${LOGROOT}/perf_summary.txt)"
  fi
fi

exit 0
