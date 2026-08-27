#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
header=$'run_accession\tlibrary_layout\tfastq_ftp\tfastq_md5'

printf '%s\n' "$header" \
  $'RUN_C\tSINGLE\thost/RUN_C.fastq.gz\tcccc' \
  $'RUN_A\tPAIRED\thost/RUN_A_1.fastq.gz;host/RUN_A_2.fastq.gz\taaaa;bbbb' \
  $'RUN_B\tSINGLE\thost/RUN_B.fastq.gz\tdddd' \
  $'RUN_A\tPAIRED\thost/RUN_A_1.fastq.gz;host/RUN_A_2.fastq.gz\taaaa;bbbb' >"$tmp/report.tsv"
"$repo/scripts/build_run_manifest.sh" SAMPLE Species Ref "$tmp/report.tsv" "$tmp/manifest.tsv"
[[ $(tail -n +2 "$tmp/manifest.tsv" | wc -l) -eq 3 ]]
[[ $(cut -f4 "$tmp/manifest.tsv" | tail -n +2 | paste -sd, -) == RUN_A,RUN_B,RUN_C ]]
awk -F '\t' 'NR>1 { if ($10 != 3 || $11 != "RUN_A,RUN_B,RUN_C") exit 1 }
  $4=="RUN_A" { if ($5!="PE" || $7==".") exit 1 }
  $4=="RUN_B" { if ($5!="SE" || $7!=".") exit 1 }' "$tmp/manifest.tsv"

printf '%s\n' "$header" \
  $'RUN_A\tSINGLE\thost/RUN_A.fastq.gz\taaaa' \
  $'RUN_A\tSINGLE\thost/RUN_A.fastq.gz\tbbbb' >"$tmp/conflict.tsv"
if "$repo/scripts/build_run_manifest.sh" S Sp Ref "$tmp/conflict.tsv" "$tmp/out" 2>/dev/null; then
  echo "conflicting run identity unexpectedly succeeded" >&2; exit 1
fi
printf '%s\n' "$header" >"$tmp/empty.tsv"
if "$repo/scripts/build_run_manifest.sh" S Sp Ref "$tmp/empty.tsv" "$tmp/out" 2>/dev/null; then
  echo "zero-run manifest unexpectedly succeeded" >&2; exit 1
fi
echo "deterministic run manifest: PASS"
