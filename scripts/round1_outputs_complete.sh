#!/usr/bin/env bash
# Authoritative Round1 published-output contract. Keep this as the single source
# used by both the Round1 Nextflow skip logic and the streaming launcher.
set -euo pipefail

[[ $# -eq 2 ]] || { echo "usage: $0 OUTDIR SAMPLE_ID" >&2; exit 2; }
outdir=${1%/}
sample=$2
root="${outdir}/${sample}/round_1"

required=(
  "numt_decoy_ref/${sample}.chrM_plus_numt.fa"
  "numt_decoy_ref/${sample}.chrM_plus_numt.fa.fai"
  "numt_decoy_ref/${sample}.original_numt.fa"
  "numt_decoy_ref/${sample}.original_numt.fa.fai"
  "numt_decoy_ref/${sample}.decoy_numt.interval_list"
  "candidate_reads/${sample}.with_mates.bam"
  "candidate_reads/${sample}.with_mates.bam.bai"
  "decoy_realign/${sample}.decoy_realign.bam"
  "decoy_realign/${sample}.decoy_realign.bam.bai"
  "numt_decoy_variant_calling/${sample}.numt_decoy.raw.vcf.gz"
  "numt_decoy_variant_calling/${sample}.numt_decoy.raw.vcf.gz.tbi"
  "numt_decoy_variant_calling/${sample}.numt_decoy.pass.split.vcf.gz"
  "numt_decoy_variant_calling/${sample}.numt_decoy.pass.split.vcf.gz.tbi"
  "consensus_numt_ref/${sample}.consensus_numt.fa"
  "consensus_numt_ref/${sample}.consensus_numt.fa.fai"
  "chrM_clean/${sample}.final_chrM.sorted.bam"
  "chrM_clean/${sample}.final_chrM.sorted.bam.bai"
  "alignment_numt_decoy/${sample}.numt_decoy.clean.cram"
  "alignment_numt_decoy/${sample}.numt_decoy.clean.cram.crai"
)

missing=0
for relative in "${required[@]}"; do
  if [[ ! -e "${root}/${relative}" ]]; then
    printf 'MISSING %s\n' "${root}/${relative}" >&2
    missing=1
  fi
done
exit "$missing"
