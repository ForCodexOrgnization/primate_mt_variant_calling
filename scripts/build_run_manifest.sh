#!/usr/bin/env bash
# Convert a resolved ENA report into the complete, deterministic run manifest
# consumed by preprocessing.nf.  No downloads are performed here.
set -euo pipefail

[[ $# -eq 5 ]] || { echo "usage: $0 SAMPLE_ID SPECIES_NAME REF_NAME RESOLVED_REPORT OUTPUT" >&2; exit 2; }
sample_id=$1 species_name=$2 ref_name=$3 report=$4 output=$5
source "$(dirname "$0")/classify_ena_fastq.sh"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
rows="$tmp/rows.tsv"
: >"$rows"

tail -n +2 "$report" | LC_ALL=C sort -t$'\t' -k1,1 -k2,2 -k3,3 -k4,4 | \
while IFS=$'\t' read -r run_id library_layout ftp_urls md5s; do
    [[ -n "$run_id" ]] || continue
    classify_ena_fastq "$run_id" "$library_layout" "$ftp_urls" "$md5s" || exit 42
    r2_url=. r2_md5=.
    if [[ "$FASTQ_LAYOUT" == PE ]]; then
        r2_url="https://$FASTQ_R2_URL"; r2_md5=$FASTQ_R2_MD5
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample_id" "$species_name" "$ref_name" "$run_id" "$FASTQ_LAYOUT" \
        "https://$FASTQ_R1_URL" "$r2_url" "$FASTQ_R1_MD5" "$r2_md5"
done >"$rows"

[[ -s "$rows" ]] || { echo "ERROR: No usable FASTQ runs found for $sample_id." >&2; exit 42; }

# Identical rows are harmless, but the same run identity must never describe
# two different inputs.  This check intentionally precedes deduplication.
awk -F '\t' '
  { key=$4; row=$0; if (key in seen && seen[key] != row) {
      printf "ERROR: DETERMINISTIC_FASTQ_FAILURE class=CONFLICTING_RUN_ID run=%s\n", key > "/dev/stderr"; exit 42
    } seen[key]=row }
' "$rows" || exit $?
LC_ALL=C sort -u -t$'\t' -k4,4 "$rows" >"$tmp/unique.tsv"
count=$(wc -l <"$tmp/unique.tsv" | tr -d ' ')
expected_ids=$(cut -f4 "$tmp/unique.tsv" | paste -sd, -)

printf 'sample_id\tspecies_name\tref_name\trun_id\tlayout\tr1_url\tr2_url\tr1_md5\tr2_md5\texpected_run_count\texpected_run_ids\n' >"$output"
awk -F '\t' -v OFS='\t' -v n="$count" -v ids="$expected_ids" '{ print $0,n,ids }' "$tmp/unique.tsv" >>"$output"
