#!/usr/bin/env bash
# Fast completion-state check for a published CRAM/CRAI pair.
set -u

cram=""
crai=""
marker=""
samtools_bin="${SAMTOOLS_BIN:-samtools}"
min_cram=1024
min_crai=16

while (($#)); do
    case "$1" in
        --cram) cram="$2"; shift 2 ;;
        --crai) crai="$2"; shift 2 ;;
        --marker) marker="$2"; shift 2 ;;
        --samtools) samtools_bin="$2"; shift 2 ;;
        --min-cram-size) min_cram="$2"; shift 2 ;;
        --min-crai-size) min_crai="$2"; shift 2 ;;
        # Deprecated compatibility options. Fast validation never uses them.
        --reference|--stability-retries|--retries|--delay|--timeout) shift 2 ;;
        *)
            echo "STATUS=UNKNOWN"
            echo "REASON=invalid_argument:$1"
            exit 2
            ;;
    esac
done

emit() {
    printf '%s\n' \
        "STATUS=$1" \
        "REASON=$2" \
        "CRAM_SIZE=${cram_size:-unknown}" \
        "CRAI_SIZE=${crai_size:-unknown}"
}

[[ -n "$cram" && -n "$crai" ]] || { emit UNKNOWN missing_required_path; exit 2; }
[[ -f "$cram" ]] || { emit INCOMPLETE cram_missing; exit 1; }
[[ -f "$crai" ]] || { emit INCOMPLETE crai_missing; exit 1; }
[[ -r "$cram" && -r "$crai" ]] || { emit UNKNOWN files_not_readable; exit 2; }

cram_size=$(stat -c%s -- "$cram" 2>/dev/null) || { emit UNKNOWN cram_stat_failed; exit 2; }
crai_size=$(stat -c%s -- "$crai" 2>/dev/null) || { emit UNKNOWN crai_stat_failed; exit 2; }

((cram_size >= min_cram)) || { emit INCOMPLETE cram_below_minimum_size; exit 1; }
((crai_size >= min_crai)) || { emit INCOMPLETE crai_below_minimum_size; exit 1; }

# Fastest path: the producer wrote this marker only after quickcheck succeeded.
if [[ -n "$marker" && -f "$marker" ]]; then
    marked_cram=$(awk -F= '$1=="cram_size" {print $2; exit}' "$marker" 2>/dev/null)
    marked_crai=$(awk -F= '$1=="crai_size" {print $2; exit}' "$marker" 2>/dev/null)
    if [[ -n "$marked_cram" && -n "$marked_crai" &&
          "$marked_cram" == "$cram_size" && "$marked_crai" == "$crai_size" ]]; then
        emit COMPLETE completion_marker_and_sizes_match
        exit 0
    fi
fi

# Legacy result: retry the lightweight quickcheck once, with no reference decode.
version_output=$("$samtools_bin" --version 2>&1)
version_rc=$?
if ((version_rc != 0)); then
    emit UNKNOWN samtools_unavailable
    exit 2
fi

quickcheck_output=""
quickcheck_rc=1

for attempt in 1 2; do
    output=$("$samtools_bin" quickcheck -v "$cram" 2>&1)
    rc=$?

    echo "QUICKCHECK_ATTEMPT=${attempt}/2 EXIT=${rc} OUTPUT=${output//$'\n'/\\n}" >&2

    if ((rc == 0)); then
        emit COMPLETE files_exist_and_quickcheck_passed
        exit 0
    fi

    quickcheck_output+=$'\n'"$output"
    quickcheck_rc=$rc
    ((attempt == 1)) && sleep 1
done

lower=${quickcheck_output,,}

if [[ "$lower" =~ truncated.file|missing.eof.block|invalid.cram|malformed|was.not.identified.as.sequence.data ]]; then
    emit INCOMPLETE confirmed_corrupt
    exit 1
fi

emit UNKNOWN quickcheck_indeterminate
exit 2
