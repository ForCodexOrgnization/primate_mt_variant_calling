#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
checker="$repo/scripts/round1_outputs_complete.sh"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
sample=ERS14600652
root="$tmp/results/$sample/round_1"
files=(
 "numt_decoy_ref/$sample.chrM_plus_numt.fa" "numt_decoy_ref/$sample.chrM_plus_numt.fa.fai"
 "numt_decoy_ref/$sample.original_numt.fa" "numt_decoy_ref/$sample.original_numt.fa.fai"
 "numt_decoy_ref/$sample.decoy_numt.interval_list"
 "candidate_reads/$sample.with_mates.bam" "candidate_reads/$sample.with_mates.bam.bai"
 "decoy_realign/$sample.decoy_realign.bam" "decoy_realign/$sample.decoy_realign.bam.bai"
 "numt_decoy_variant_calling/$sample.numt_decoy.raw.vcf.gz" "numt_decoy_variant_calling/$sample.numt_decoy.raw.vcf.gz.tbi"
 "numt_decoy_variant_calling/$sample.numt_decoy.pass.split.vcf.gz" "numt_decoy_variant_calling/$sample.numt_decoy.pass.split.vcf.gz.tbi"
 "consensus_numt_ref/$sample.consensus_numt.fa" "consensus_numt_ref/$sample.consensus_numt.fa.fai"
 "chrM_clean/$sample.final_chrM.sorted.bam" "chrM_clean/$sample.final_chrM.sorted.bam.bai"
 "alignment_numt_decoy/$sample.numt_decoy.clean.cram" "alignment_numt_decoy/$sample.numt_decoy.clean.cram.crai"
)
for f in "${files[@]}"; do mkdir -p "$root/${f%/*}"; : >"$root/$f"; done
vcf_tree="$root/mtdna_variant_calling/cromwell-executions/MitochondriaMultiSamplePipeline/fake-workflow-id/call-SplitMultiAllelicSites/shard-0/execution"
mkdir -p "$vcf_tree"
printf '##fileformat=VCFv4.2\n' >"$vcf_tree/$sample.numt_decoy.clean.final.split.vcf"

# Production regression: the obsolete launcher-only WDL VCF is absent, but the
# exact current workflow skip contract is complete.
[[ ! -e "$tmp/results/$sample/round_1_variant_calling_decoy/$sample.numt_decoy.clean.final.split.vcf" ]]
"$checker" "$tmp/results" "$sample"

# Mirror FIND_ROUND1_OUTPUTS' exact-name recursive discovery and prove the
# current Cromwell tree is selected without creating the legacy directory.
mapfile -t round2_hits < <(find "$root/mtdna_variant_calling" -type f -name "$sample.numt_decoy.clean.final.split.vcf" | sort)
[[ ${#round2_hits[@]} -eq 1 ]]
[[ "${round2_hits[0]}" == "$vcf_tree/$sample.numt_decoy.clean.final.split.vcf" ]]

# A genuinely required current index makes the sample incomplete.
rm "$root/candidate_reads/$sample.with_mates.bam.bai"
if "$checker" "$tmp/results" "$sample" 2>/dev/null; then
  echo 'missing required Round1 BAM index was accepted' >&2; exit 1
fi

echo 'round1 authoritative completion tests: PASS'
