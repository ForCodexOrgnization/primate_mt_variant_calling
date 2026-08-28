#!/usr/bin/env bash
# Retry-safe, serialized update of the global successful-sample provenance index.
set -euo pipefail

[[ $# -eq 5 ]] || { echo "usage: $0 SUMMARY_TSV SAMPLE_ID SPECIES_NAME REF_NAME RUN_IDS" >&2; exit 2; }
summary=$1 sample_id=$2 species_name=$3 ref_name=$4 run_ids=$5
mkdir -p "$(dirname "$summary")"

normalized_ids=$(printf '%s\n' "$run_ids" | tr ',' '\n' | sed '/^$/d' | LC_ALL=C sort -u | paste -sd, -)
[[ -n "$normalized_ids" ]] || { echo "ERROR: empty used run set for $sample_id" >&2; exit 2; }
count=$(awk -F, '{print NF}' <<<"$normalized_ids")

exec 9>"${summary}.lock"
flock 9
tmp=$(mktemp "$(dirname "$summary")/.used_short_read_runs.tsv.XXXXXX")
trap 'rm -f "$tmp"' EXIT
printf 'sample_id\tspecies_name\tref_name\tshort_read_run_count\tshort_read_run_ids\n' >"$tmp"
if [[ -f "$summary" ]]; then
    awk -F '\t' -v sample="$sample_id" 'NR > 1 && $1 != sample' "$summary" >>"$tmp"
fi
printf '%s\t%s\t%s\t%s\t%s\n' "$sample_id" "$species_name" "$ref_name" "$count" "$normalized_ids" >>"$tmp"
{ head -n1 "$tmp"; tail -n +2 "$tmp" | LC_ALL=C sort -t$'\t' -k1,1; } >"${tmp}.sorted"
mv "${tmp}.sorted" "$summary"
