#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
helper="$repo/scripts/sam_platform.sh"
[[ $($helper map ILLUMINA) == ILLUMINA ]]
[[ $($helper map DNBSEQ) == DNBSEQ ]]
[[ $($helper map BGISEQ) == DNBSEQ ]]
if "$helper" map PACBIO_SMRT >/dev/null 2>&1; then exit 1; fi

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir "$tmp/bin"
cat >"$tmp/bin/samtools" <<'EOF'
#!/usr/bin/env bash
case "$3" in
  illumina) printf '@HD\tVN:1.6\n@RG\tID:a\tPL:ILLUMINA\n@RG\tID:b\tPL:ILLUMINA\n' ;;
  dnbseq) printf '@RG\tID:a\tPL:DNBSEQ\n' ;;
  mixed) printf '@RG\tID:a\tPL:DNBSEQ\n@RG\tID:b\tPL:ILLUMINA\n' ;;
  empty) printf '@HD\tVN:1.6\n@RG\tID:a\n' ;;
  bad) printf '@RG\tID:a\tPL:PACBIO\n' ;;
esac
EOF
chmod +x "$tmp/bin/samtools"
PATH="$tmp/bin:$PATH"
[[ $($helper derive illumina) == ILLUMINA ]]
[[ $($helper derive dnbseq) == DNBSEQ ]]
for input in mixed empty bad; do if "$helper" derive "$input" >/dev/null 2>&1; then exit 1; fi; done

# Source contracts cover all preprocessing alignment routes and both NUMT rounds.
! rg -q 'PL:ILLUMINA' "$repo/preprocessing.nf" "$repo/primate_pipeline_numt_decoy_round1.nf" "$repo/primate_pipeline_round2_consensus_NUMT.nf"
[[ $(rg -c 'PL:\$\{meta.sam_platform\}' "$repo/preprocessing.nf") -eq 3 ]]
rg -q 'sam_platform.sh.*derive.*bam_with_mates' "$repo/primate_pipeline_numt_decoy_round1.nf"
rg -q 'sam_platform.sh.*derive.*round1_bam' "$repo/primate_pipeline_round2_consensus_NUMT.nf"
rg -q 'instrument_platform.*row.instrument_platform' "$repo/preprocessing.nf"
# The sole remaining literal is explicitly scoped to the legacy Illumina-only workflow.
[[ $(rg -l 'PL:ILLUMINA' "$repo" --glob '*.nf' | wc -l) -eq 1 ]]
rg -q 'Legacy SRA-only workflow' "$repo/primate_pipeline.nf"
# Deterministic markers no longer require run= and reason= suffixes.
rg -q 'DETERMINISTIC_FASTQ_FAILURE.*class=' "$repo/launch_pipeline_streaming_per_sample.sh"
! rg -q 'class=.*run=.*reason=' "$repo/launch_pipeline_streaming_per_sample.sh"
echo "SAM platform mapping and propagation: PASS"
