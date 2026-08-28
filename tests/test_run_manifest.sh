#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
header=$'run_accession\tlibrary_layout\tfastq_ftp\tfastq_md5\tinstrument_platform\tinstrument_model'

pe() { printf '%s\t%s\th/%s_1.fastq.gz;h/%s_2.fastq.gz\t%s1;%s2\t%s\t%s\n' "$1" "${2:-PAIRED}" "$1" "$1" "$1" "$1" "$3" "$4"; }
se() { printf '%s\t%s\th/%s.fastq.gz\t%s\t%s\t%s\n' "$1" "$2" "$1" "$1" "$3" "$4"; }
run_ok() {
    local name=$1; shift
    { printf '%s\n' "$header"; printf '%s\n' "$@"; } >"$tmp/$name.report"
    "$repo/scripts/build_run_manifest.sh" "$name" Species Ref "$tmp/$name.report" "$tmp/$name.manifest" "$tmp/$name.excluded"
}
run_fail() {
    local name=$1; shift
    { printf '%s\n' "$header"; printf '%s\n' "$@"; } >"$tmp/$name.report"
    set +e
    "$repo/scripts/build_run_manifest.sh" "$name" Species Ref "$tmp/$name.report" "$tmp/$name.manifest" "$tmp/$name.excluded" 2>"$tmp/$name.err"
    local rc=$?
    set -e
    [[ $rc -eq 42 ]]
    grep -q 'DETERMINISTIC_FASTQ_FAILURE class=NO_SUPPORTED_SHORT_READ_RUNS' "$tmp/$name.err"
}

# 1: Illumina PE only; all eligible runs are deterministic and complete.
run_ok ILLUMINA_ONLY "$(pe I2 PAIRED ILLUMINA NovaSeq)" "$(pe I1 SINGLE ILLUMINA HiSeq)"
[[ $(tail -n +2 "$tmp/ILLUMINA_ONLY.manifest" | wc -l) -eq 2 ]]
awk -F '\t' 'NR>1 {if ($5!="ILLUMINA" || $7!="PE" || $12!=2 || $13!="I1,I2") exit 1}' "$tmp/ILLUMINA_ONLY.manifest"

# 2-3: Actual layout, not metadata layout, governs eligibility.
run_ok ILLUMINA_MIX "$(pe IPE PAIRED ILLUMINA NovaSeq)" "$(se ISE PAIRED ILLUMINA NovaSeq)"
[[ $(tail -n +2 "$tmp/ILLUMINA_MIX.manifest" | cut -f4) == IPE ]]
grep -q $'ISE\tILLUMINA\tNovaSeq\tUNSUPPORTED_FASTQ_LAYOUT_SE' "$tmp/ILLUMINA_MIX.excluded"
run_fail ILLUMINA_SE "$(se ISE SINGLE ILLUMINA HiSeq)"

# 4-6: DNBSEQ and BGISEQ are supported, but only with actual paired files.
run_ok DNB_ONLY "$(pe D1 PAIRED DNBSEQ DNBSEQ-G400)"
awk -F '\t' 'NR==2 {exit !($5=="DNBSEQ" && $7=="PE")}' "$tmp/DNB_ONLY.manifest"
run_ok BGI_ONLY "$(pe B1 PAIRED BGISEQ BGISEQ-500)"
awk -F '\t' 'NR==2 {exit !($5=="BGISEQ" && $7=="PE")}' "$tmp/BGI_ONLY.manifest"
run_fail BGI_SE "$(se BSE PAIRED BGISEQ BGISEQ-500)"
grep -q 'UNSUPPORTED_FASTQ_LAYOUT_SE' "$tmp/BGI_SE.excluded"

# 7: Illumina wins over eligible BGISEQ; unsupported PacBio is never classified.
pacbio=$'PB\tSINGLE\th/PB_subreads.fastq.gz\tpb\tPACBIO_SMRT\tSequel II'
run_ok ALL_MIX "$(pe I1 PAIRED ILLUMINA NovaSeq)" "$(pe B1 PAIRED BGISEQ BGISEQ-500)" "$pacbio"
[[ $(tail -n +2 "$tmp/ALL_MIX.manifest" | cut -f4) == I1 ]]
grep -q $'B1\tBGISEQ\tBGISEQ-500\tLOWER_PRIORITY_SUPPORTED_PLATFORM' "$tmp/ALL_MIX.excluded"
grep -q $'PB\tPACBIO_SMRT\tSequel II\tUNSUPPORTED_SEQUENCING_PLATFORM' "$tmp/ALL_MIX.excluded"

# 8: Audit-shaped sample: Illumina PE remains; BGISEQ SE and PacBio do not.
run_ok AUDIT_SAMPLE "$(pe I1 PAIRED ILLUMINA NovaSeq)" "$(se BSE PAIRED BGISEQ BGISEQ-500)" "$pacbio"
[[ $(tail -n +2 "$tmp/AUDIT_SAMPLE.manifest" | cut -f4) == I1 ]]
grep -q $'BSE\tBGISEQ\tBGISEQ-500\tUNSUPPORTED_FASTQ_LAYOUT_SE' "$tmp/AUDIT_SAMPLE.excluded"

# 9-10: DNBSEQ wins over BGISEQ; BGISEQ wins when it is the only eligible PE platform.
run_ok DNB_BGI "$(pe B1 PAIRED BGISEQ BGI)" "$(pe D1 PAIRED DNBSEQ DNB)"
[[ $(tail -n +2 "$tmp/DNB_BGI.manifest" | cut -f5 | sort -u) == DNBSEQ ]]
grep -q 'LOWER_PRIORITY_SUPPORTED_PLATFORM' "$tmp/DNB_BGI.excluded"
run_ok BGI_PB "$(pe B1 PAIRED BGISEQ BGI)" "$pacbio"
[[ $(tail -n +2 "$tmp/BGI_PB.manifest" | cut -f5) == BGISEQ ]]

# 11-13: Unsupported-only samples terminate deterministically without classification.
run_fail PACBIO_ONLY "$pacbio"
run_fail ONT_ONLY $'ONT\tSINGLE\th/ONT.fastq.gz\to\tOXFORD_NANOPORE\tMinION'
run_fail LS454_ONLY $'LS\tSINGLE\th/LS.fastq.gz\tl\tLS454\t454 GS FLX'
for n in PACBIO_ONLY ONT_ONLY LS454_ONLY; do grep -q 'UNSUPPORTED_SEQUENCING_PLATFORM' "$tmp/$n.excluded"; done

# 14: All selected-platform runs are retained, exact, and platform homogeneous.
run_ok MULTI "$(pe D2 PAIRED DNBSEQ DNB)" "$(pe B1 PAIRED BGISEQ BGI)" "$(pe D1 PAIRED DNBSEQ DNB)"
awk -F '\t' 'NR>1 {if ($5!="DNBSEQ" || $12!=2 || $13!="D1,D2") exit 1}' "$tmp/MULTI.manifest"
[[ $(tail -n +2 "$tmp/MULTI.manifest" | cut -f5 | sort -u | wc -l) -eq 1 ]]

# 15: Conflicting duplicate eligible run identity keeps the public failure.
{ printf '%s\n' "$header"; pe DUP PAIRED ILLUMINA HiSeq; printf '%s\n' $'DUP\tPAIRED\th/DUP_1.fastq.gz;h/DUP_2.fastq.gz\tother1;other2\tILLUMINA\tHiSeq'; } >"$tmp/conflict.report"
set +e
"$repo/scripts/build_run_manifest.sh" CONFLICT Species Ref "$tmp/conflict.report" "$tmp/conflict.manifest" 2>"$tmp/conflict.err"
rc=$?
set -e
[[ $rc -eq 42 ]]
grep -q 'class=CONFLICTING_RUN_ID run=DUP' "$tmp/conflict.err"

echo "deterministic run manifest platform policy: PASS"
