#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
Usage:
  aggregate_kmc_stats_and_totals_1.sh \
    --metadata <metadata.tsv> \
    --outdir <results_dir> \
    [--outfile-samples <path>] \
    [--outfile-pooled  <path>] \
    [--logfile <path>]
USAGE
  exit 1
}

META="metadata.tsv"
OUTDIR="results"
OUT_SAMPLES=""
OUT_POOLED=""
LOGFILE=""

BIN_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
ROOT_DIR="$(cd -- "${BIN_DIR}/.." >/dev/null 2>&1 && pwd)"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --metadata) META="$2"; shift 2;;
    --outdir)   OUTDIR="$2"; shift 2;;
    --outfile-samples) OUT_SAMPLES="$2"; shift 2;;
    --outfile-pooled)  OUT_POOLED="$2"; shift 2;;
    --logfile)  LOGFILE="$2"; shift 2;;
    -h|--help)  usage;;
    *) echo "Unknown arg: $1" >&2; usage;;
  esac
done

[[ -f "$META" ]] || { echo "[aggregate] ERR: metadata not found: $META"; exit 1; }
[[ -n "$LOGFILE" ]] || LOGFILE="${ROOT_DIR}/logging/build_pooled_kmc_files.log"
[[ -f "$LOGFILE" ]] || { echo "[aggregate] ERR: log not found: $LOGFILE"; exit 1; }

STATDIR="${OUTDIR%/}/statistics"
mkdir -p "${STATDIR}"
[[ -n "$OUT_SAMPLES" ]] || OUT_SAMPLES="${STATDIR}/kmer_stats_per_sample.tsv"
[[ -n "$OUT_POOLED"  ]] || OUT_POOLED="${STATDIR}/kmer_totals.tsv"

head1="$(head -n1 "$META")"
if   grep -q $'\t' <<<"$head1"; then DELIM=$'\t'
elif grep -q ';'   <<<"$head1"; then DELIM=';'
elif grep -q ','   <<<"$head1"; then DELIM=','
else                                 DELIM='[ \t]+'
fi

TMP_LIST="$(mktemp)"
awk -v FS="${DELIM}" '
/^[ \t\r]*#/ || /^[ \t\r]*$/ { next }
NR==1{
  for(i=1;i<=NF;i++){ c=$i; gsub(/^[ \t\r]+|[ \t\r]+$/,"",c); lc=tolower(c); h[lc]=i }
  if(!("trait" in h)){ print "[aggregate] ERR: metadata needs column trait" > "/dev/stderr"; exit 1 }
  if(!("sample_id" in h)){ print "[aggregate] ERR: metadata needs column sample_id" > "/dev/stderr"; exit 1 }
  next
}
{
  sid=$h["sample_id"]; tr=$h["trait"]
  gsub(/^[ \t\r]+|[ \t\r]+$/,"",sid); gsub(/^[ \t\r]+|[ \t\r]+$/,"",tr)
  if(sid!="" && tr!="") print sid "\t" tr
}' "$META" | sort -u > "${TMP_LIST}"

echo -e "sample_id\ttrait\tbelow_min\tabove_max\tunique_kmers\tunique_counted_kmers\ttotal_kmers\ttotal_reads\ttotal_super_kmers" > "${OUT_SAMPLES}"

extract_sample_stats() {
  sid="$1"
  awk -v SID="$sid" '
    BEGIN{inblk=0; stats=0; seen=0; below=above=uniqa=uniqc=total=reads=super=""}
    /^\[create\] /{ if ($0 ~ ("--sample[[:space:]]*" SID "([[:space:]]|$)")) { inblk=1; next } }
    inblk && /^Stats:/ { stats=1; next }
    inblk && stats && /No\. of k-mers below min\. threshold/ { s=$0; sub(/.*:/,"",s); gsub(/[[:space:]]/,"",s); below=s }
    inblk && stats && /No\. of k-mers above max\. threshold/ { s=$0; sub(/.*:/,"",s); gsub(/[[:space:]]/,"",s); above=s }
    inblk && stats && /No\. of unique k-mers[[:space:]]*:/   { s=$0; sub(/.*:/,"",s); gsub(/[[:space:]]/,"",s); uniqa=s }
    inblk && stats && /No\. of unique counted k-mers/        { s=$0; sub(/.*:/,"",s); gsub(/[[:space:]]/,"",s); uniqc=s }
    inblk && stats && /Total no\. of k-mers[[:space:]]*:/    { s=$0; sub(/.*:/,"",s); gsub(/[[:space:]]/,"",s); total=s }
    inblk && stats && /Total no\. of reads[[:space:]]*:/     { s=$0; sub(/.*:/,"",s); gsub(/[[:space:]]/,"",s); reads=s }
    inblk && stats && /Total no\. of super-k-mers[[:space:]]*:/ { s=$0; sub(/.*:/,"",s); gsub(/[[:space:]]/,"",s); super=s }
    /^\[create\] DONE[[:space:]]*/{ if (inblk){ inblk=0; stats=0; seen=1 } }
    END{
      if(seen){
        if(below=="") below=0; if(above=="") above=0; if(uniqa=="") uniqa=0;
        if(uniqc=="") uniqc=0; if(total=="") total=0; if(reads=="") reads=0; if(super=="") super=0;
        print below "\t" above "\t" uniqa "\t" uniqc "\t" total "\t" reads "\t" super;
      }
    }
  ' "${LOGFILE}"
}

TMP_TOT="$(mktemp)"
: > "${TMP_TOT}"

while IFS=$'\t' read -r SID TRAIT; do
  [[ -z "$SID" ]] && continue
  stats_line="$(extract_sample_stats "$SID" || true)"
  if [[ -z "$stats_line" ]]; then
    echo -e "${SID}\t${TRAIT}\t0\t0\t0\t0\t0\t0\t0" >> "${OUT_SAMPLES}"
    echo -e "${TRAIT}\t0\t0\t0" >> "${TMP_TOT}"
    continue
  fi
  IFS=$'\t' read -r below above uniqa uniqc total reads super <<EOF
${stats_line}
EOF
  echo -e "${SID}\t${TRAIT}\t${below}\t${above}\t${uniqa}\t${uniqc}\t${total}\t${reads}\t${super}" >> "${OUT_SAMPLES}"
  echo -e "${TRAIT}\t${total}\t${reads}\t${super}" >> "${TMP_TOT}"
done < "${TMP_LIST}"

rm -f "${TMP_LIST}"

{
  echo -e "trait\ttotal_kmers_sum\ttotal_reads_sum\ttotal_super_kmers_sum"
  awk -F'\t' '{ a[$1,1]+=$2; a[$1,2]+=$3; a[$1,3]+=$4; k[$1]=1 } END{ for(t in k) printf "%s\t%d\t%d\t%d\n", t, a[t,1]+0, a[t,2]+0, a[t,3]+0 }' "${TMP_TOT}" | sort
} > "${OUT_POOLED}"

rm -f "${TMP_TOT}"
echo "[aggregate] Wrote ${OUT_SAMPLES}"
echo "[aggregate] Wrote ${OUT_POOLED}"
