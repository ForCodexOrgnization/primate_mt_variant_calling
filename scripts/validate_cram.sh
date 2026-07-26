#!/usr/bin/env bash
# Validate a published CRAM without confusing validator/storage failures with corruption.
set -u

cram= crai= reference= marker= samtools_bin="${SAMTOOLS_BIN:-samtools}"
retries=3 stability_retries=3 delay=10 timeout_seconds=120 min_cram=1024 min_crai=16
while (($#)); do
  case "$1" in
    --cram) cram=$2; shift 2;; --crai) crai=$2; shift 2;;
    --reference) reference=$2; shift 2;; --marker) marker=$2; shift 2;;
    --samtools) samtools_bin=$2; shift 2;; --retries) retries=$2; shift 2;;
    --delay) delay=$2; shift 2;; --timeout) timeout_seconds=$2; shift 2;;
    --stability-retries) stability_retries=$2; shift 2;;
    --min-cram-size) min_cram=$2; shift 2;; --min-crai-size) min_crai=$2; shift 2;;
    *) echo "STATUS=UNKNOWN"; echo "REASON=invalid_argument:$1"; exit 2;;
  esac
done

emit() { printf '%s\n' "STATUS=$1" "REASON=$2" "CRAM_SIZE=${cram_size:-unknown}" "CRAI_SIZE=${crai_size:-unknown}" "SAMTOOLS_VERSION=${samtools_version:-unknown}"; }
unknown_re='input/output error|stale file handle|resource temporarily unavailable|permission denied|transport endpoint|timed out|timeout|no such device|temporarily unavailable'
corrupt_re='was not identified as sequence data|truncated file|missing eof block|malformed|invalid cram|error reading file|cram.*corrupt'

[[ -n "$cram" && -n "$crai" ]] || { emit UNKNOWN missing_required_path; exit 2; }
[[ -e "$cram" ]] || { emit INCOMPLETE cram_missing; exit 1; }
[[ -e "$crai" ]] || { emit INCOMPLETE crai_missing; exit 1; }
[[ -f "$cram" && -r "$cram" && -f "$crai" && -r "$crai" ]] || { emit UNKNOWN files_not_readable; exit 2; }

stable=false
for ((attempt=1; attempt<=stability_retries; attempt++)); do
  size1=$(stat -c%s -- "$cram" 2>&1) || { emit UNKNOWN "cram_stat_failed:${size1//$'\n'/ }"; exit 2; }
  isize1=$(stat -c%s -- "$crai" 2>&1) || { emit UNKNOWN "crai_stat_failed:${isize1//$'\n'/ }"; exit 2; }
  sleep "$delay"
  size2=$(stat -c%s -- "$cram" 2>&1) || { emit UNKNOWN "cram_stat_failed:${size2//$'\n'/ }"; exit 2; }
  isize2=$(stat -c%s -- "$crai" 2>&1) || { emit UNKNOWN "crai_stat_failed:${isize2//$'\n'/ }"; exit 2; }
  cram_size=$size2; crai_size=$isize2
  echo "STABILITY_ATTEMPT=$attempt/$stability_retries CRAM=$size1->$size2 CRAI=$isize1->$isize2" >&2
  if [[ $size1 == "$size2" && $isize1 == "$isize2" ]]; then
    if ((size2 < min_cram)); then emit INCOMPLETE cram_below_minimum_size; exit 1; fi
    if ((isize2 < min_crai)); then emit INCOMPLETE crai_below_minimum_size; exit 1; fi
    stable=true; break
  fi
done
[[ $stable == true ]] || { emit UNKNOWN files_not_stable; exit 2; }

version_output=$("$samtools_bin" --version 2>&1); version_rc=$?
((version_rc == 0)) || { emit UNKNOWN "samtools_unavailable:${version_output//$'\n'/ }"; exit 2; }
samtools_version=$(printf '%s' "$version_output" | head -n1)

# A marker is an optimization/safety record, never a substitute for validation.
if [[ -n "$marker" && -f "$marker" ]]; then
  marked_cram=$(awk -F= '$1=="cram_size" {print $2; exit}' "$marker" 2>/dev/null)
  marked_crai=$(awk -F= '$1=="crai_size" {print $2; exit}' "$marker" 2>/dev/null)
  if [[ -z "$marked_cram" || -z "$marked_crai" || "$marked_cram" != "$cram_size" || "$marked_crai" != "$crai_size" ]]; then
    emit UNKNOWN marker_size_mismatch; exit 2
  fi
fi

all_output= all_corrupt=true
for ((attempt=1; attempt<=retries; attempt++)); do
  output=$(timeout --signal=KILL "$timeout_seconds" "$samtools_bin" quickcheck -v "$cram" 2>&1); rc=$?
  printf 'QUICKCHECK_ATTEMPT=%s/%s EXIT=%s\nCOMMAND=%q quickcheck -v %q\nOUTPUT=%s\nCRAM_SIZE=%s CRAI_SIZE=%s\n' \
    "$attempt" "$retries" "$rc" "$samtools_bin" "$cram" "$output" "$cram_size" "$crai_size" >&2
  if ((rc == 0)); then all_corrupt=false; quickcheck_attempts=$attempt; break; fi
  all_output+=$'\n'$output
  lower=${output,,}
  [[ $lower =~ $corrupt_re ]] || all_corrupt=false
  ((rc == 124 || rc == 137)) && all_corrupt=false
  ((attempt < retries)) && sleep "$delay"
done
if ((rc != 0)); then
  lower=${all_output,,}
  if [[ $all_corrupt == true ]]; then emit INCOMPLETE confirmed_corrupt_after_retries; exit 1; fi
  [[ $lower =~ $unknown_re ]] && reason=transient_or_storage_error || reason=quickcheck_indeterminate
  emit UNKNOWN "$reason"; exit 2
fi

idx=("$samtools_bin" idxstats)
if [[ -n "$reference" ]]; then
  [[ -r "$reference" && -r "${reference}.fai" ]] || { emit UNKNOWN reference_or_fai_unavailable; exit 2; }
  idx+=(
    --input-fmt-option
    "reference=${reference}"
  )
fi
idx+=("$cram")
idx_output=$(timeout --signal=KILL "$timeout_seconds" "${idx[@]}" 2>&1); idx_rc=$?
echo "IDXSTATS_EXIT=$idx_rc COMMAND=${idx[*]} OUTPUT=${idx_output//$'\n'/\\n}" >&2
if ((idx_rc != 0)); then
  if [[ $idx_output == *"Usage: samtools idxstats"* ]]; then
    emit UNKNOWN idxstats_command_invalid
  else
    emit UNKNOWN idxstats_failed
  fi
  exit 2
fi
emit COMPLETE quickcheck_and_idxstats_passed
echo "QUICKCHECK_ATTEMPTS=$quickcheck_attempts"
exit 0
