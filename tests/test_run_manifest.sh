#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
header=$'run_accession\tlibrary_layout\tfastq_ftp\tfastq_md5\tinstrument_platform\tinstrument_model'

printf '%s\n' "$header" \
  $'RUN_C\tSINGLE\thost/RUN_C.fastq.gz\tcccc\tILLUMINA\tNovaSeq 6000' \
  $'RUN_A\tPAIRED\thost/RUN_A_1.fastq.gz;host/RUN_A_2.fastq.gz\taaaa;bbbb\tILLUMINA\tHiSeq 2000' \
  $'RUN_B\tSINGLE\thost/RUN_B.fastq.gz\tdddd\tILLUMINA\tGenome Analyzer' \
  $'RUN_A\tPAIRED\thost/RUN_A_1.fastq.gz;host/RUN_A_2.fastq.gz\taaaa;bbbb\tILLUMINA\tHiSeq 2000' >"$tmp/report.tsv"
"$repo/scripts/build_run_manifest.sh" SAMPLE Species Ref "$tmp/report.tsv" "$tmp/manifest.tsv"
[[ $(tail -n +2 "$tmp/manifest.tsv" | wc -l) -eq 3 ]]
[[ $(cut -f4 "$tmp/manifest.tsv" | tail -n +2 | paste -sd, -) == RUN_A,RUN_B,RUN_C ]]
awk -F '\t' 'NR>1 { if ($12 != 3 || $13 != "RUN_A,RUN_B,RUN_C") exit 1 }
  $4=="RUN_A" { if ($7!="PE" || $9==".") exit 1 }
  $4=="RUN_B" { if ($7!="SE" || $9!=".") exit 1 }' "$tmp/manifest.tsv"

printf '%s\n' "$header" \
  $'RUN_A\tSINGLE\thost/RUN_A.fastq.gz\taaaa\tILLUMINA\tHiSeq' \
  $'RUN_A\tSINGLE\thost/RUN_A.fastq.gz\tbbbb\tILLUMINA\tHiSeq' >"$tmp/conflict.tsv"
if "$repo/scripts/build_run_manifest.sh" S Sp Ref "$tmp/conflict.tsv" "$tmp/out" 2>/dev/null; then
  echo "conflicting run identity unexpectedly succeeded" >&2; exit 1
fi
printf '%s\n' "$header" >"$tmp/empty.tsv"
if "$repo/scripts/build_run_manifest.sh" S Sp Ref "$tmp/empty.tsv" "$tmp/out" 2>/dev/null; then
  echo "zero-run manifest unexpectedly succeeded" >&2; exit 1
fi

# Unsupported platforms are filtered before filename/layout classification.
printf '%s\n' "$header" \
  $'SRR17126582\tSINGLE\th/SRR17126582_subreads.fastq.gz\tdeadbeef\tPACBIO_SMRT\tPacBio RS II' >"$tmp/pacbio.tsv"
set +e
"$repo/scripts/build_run_manifest.sh" SRS11219889 Sp Ref "$tmp/pacbio.tsv" "$tmp/out" "$tmp/pacbio.excluded.tsv" 2>"$tmp/pacbio.err"
rc=$?
set -e
[[ $rc -eq 42 ]]
grep -q 'class=NO_SUPPORTED_SHORT_READ_RUNS sample=SRS11219889' "$tmp/pacbio.err"
! grep -q 'AMBIGUOUS_FASTQ_LAYOUT' "$tmp/pacbio.err"
grep -q $'SRS11219889\tSRR17126582\tPACBIO_SMRT\tPacBio RS II\tUNSUPPORTED_SEQUENCING_PLATFORM' "$tmp/pacbio.excluded.tsv"

printf '%s\n' "$header" \
  $'SRR2\tSINGLE\th/SRR2.fastq.gz\tm2\tILLUMINA\tNovaSeq 6000' \
  $'SRRPB\tSINGLE\th/SRRPB_subreads.fastq.gz\tmp\tPACBIO_SMRT\tPacBio RS II' \
  $'SRR1\tPAIRED\th/SRR1_1.fastq.gz;h/SRR1_2.fastq.gz\tm1a;m1b\tILLUMINA\tHiSeq 2000' >"$tmp/mixed.tsv"
"$repo/scripts/build_run_manifest.sh" MIXED Sp Ref "$tmp/mixed.tsv" "$tmp/mixed.manifest.tsv" "$tmp/mixed.excluded.tsv"
[[ $(tail -n +2 "$tmp/mixed.manifest.tsv" | wc -l) -eq 2 ]]
awk -F '\t' 'NR>1 { if ($5 != "ILLUMINA" || $12 != 2 || $13 != "SRR1,SRR2") exit 1 }' "$tmp/mixed.manifest.tsv"
grep -q $'MIXED\tSRRPB\tPACBIO_SMRT\tPacBio RS II\tUNSUPPORTED_SEQUENCING_PLATFORM' "$tmp/mixed.excluded.tsv"
echo "deterministic run manifest: PASS"
