#!/usr/bin/env bash
# Convert a resolved ENA report into the complete, deterministic run manifest
# consumed by preprocessing.nf.  No downloads are performed here.
set -euo pipefail

[[ $# -ge 5 && $# -le 6 ]] || { echo "usage: $0 SAMPLE_ID SPECIES_NAME REF_NAME RESOLVED_REPORT OUTPUT [EXCLUDED_RUNS_OUTPUT]" >&2; exit 2; }
sample_id=$1 species_name=$2 ref_name=$3 report=$4 output=$5
excluded_output=${6:-"${output%/*}/${sample_id}.excluded_runs.tsv"}
source "$(dirname "$0")/classify_ena_fastq.sh"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
candidates="$tmp/candidates.tsv"
: >"$candidates"
printf 'sample_id\trun_id\tinstrument_platform\tinstrument_model\texclusion_reason\n' >"$excluded_output"

# ENA library_layout is deliberately not consulted when deciding eligibility:
# classify_ena_fastq inspects the actual downloadable file set.
while IFS=$'\t' read -r run_id library_layout ftp_urls md5s instrument_platform instrument_model; do
    [[ -n "$run_id" ]] || continue
    case "$instrument_platform" in
        ILLUMINA|DNBSEQ|BGISEQ) ;;
        *)
            printf '%s\t%s\t%s\t%s\t%s\n' "$sample_id" "$run_id" "$instrument_platform" "$instrument_model" \
                UNSUPPORTED_SEQUENCING_PLATFORM >>"$excluded_output"
            continue
            ;;
    esac

    classify_ena_fastq "$run_id" "$library_layout" "$ftp_urls" "$md5s" || exit 42
    if [[ "$FASTQ_LAYOUT" != PE ]]; then
        printf '%s\t%s\t%s\t%s\t%s\n' "$sample_id" "$run_id" "$instrument_platform" "$instrument_model" \
            UNSUPPORTED_FASTQ_LAYOUT_SE >>"$excluded_output"
        continue
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample_id" "$species_name" "$ref_name" "$run_id" "$instrument_platform" "$instrument_model" PE \
        "https://$FASTQ_R1_URL" "https://$FASTQ_R2_URL" "$FASTQ_R1_MD5" "$FASTQ_R2_MD5" >>"$candidates"
done < <(tail -n +2 "$report" | LC_ALL=C sort -t$'\t' -k1,1 -k2,2 -k3,3 -k4,4 -k5,5 -k6,6)

# Identical rows are harmless, but the same run identity must never describe
# two different eligible inputs. This intentionally precedes deduplication and
# sample-level platform selection.
awk -F '\t' '
  { key=$4; row=$0; if (key in seen && seen[key] != row) {
      printf "ERROR: DETERMINISTIC_FASTQ_FAILURE class=CONFLICTING_RUN_ID run=%s\n", key > "/dev/stderr"; exit 42
    } seen[key]=row }
' "$candidates" || exit $?
LC_ALL=C sort -u -t$'\t' -k4,4 "$candidates" >"$tmp/unique.tsv"

selected_platform=
for platform in ILLUMINA DNBSEQ BGISEQ; do
    if awk -F '\t' -v p="$platform" '$5 == p { found=1; exit } END { exit !found }' "$tmp/unique.tsv"; then
        selected_platform=$platform
        break
    fi
done

if [[ -z "$selected_platform" ]]; then
    excluded=$(tail -n +2 "$excluded_output" | awk -F '\t' '{printf "%s%s/%s/%s", (NR>1 ? "," : ""), $2, $3, $4}')
    printf 'ERROR: DETERMINISTIC_FASTQ_FAILURE class=NO_SUPPORTED_SHORT_READ_RUNS sample=%s excluded_runs=%s\n' "$sample_id" "${excluded:-none}" >&2
    exit 42
fi

awk -F '\t' -v OFS='\t' -v selected="$selected_platform" -v excluded="$excluded_output" '
  $5 == selected { print; next }
  { print $1,$4,$5,$6,"LOWER_PRIORITY_SUPPORTED_PLATFORM" >> excluded }
' "$tmp/unique.tsv" >"$tmp/selected.tsv"

count=$(wc -l <"$tmp/selected.tsv" | tr -d ' ')
expected_ids=$(cut -f4 "$tmp/selected.tsv" | paste -sd, -)
printf 'sample_id\tspecies_name\tref_name\trun_id\tinstrument_platform\tinstrument_model\tlayout\tr1_url\tr2_url\tr1_md5\tr2_md5\texpected_run_count\texpected_run_ids\n' >"$output"
awk -F '\t' -v OFS='\t' -v n="$count" -v ids="$expected_ids" '{ print $0,n,ids }' "$tmp/selected.tsv" >>"$output"
