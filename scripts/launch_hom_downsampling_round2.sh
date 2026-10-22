#!/usr/bin/env bash
# Run authoritative Round2 on one completed Round1 HOM downsampling pilot.
# bash scripts/launch_hom_downsampling_round2.sh ERS12091908 100 101
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CALL_REPO="$(cd "$SCRIPT_DIR/.." && pwd)"
ANALYSIS_ROOT="${HOM_ANALYSIS_ROOT:-/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only}"
BASE="$ANALYSIS_ROOT/results/qc_coverage_only/downsampling_pilot"
SAMPLE="${1:?Provide sample ID}"; TARGET="${2:?Provide raw depth target}"; SEED="${3:?Provide seed}"
[[ "$SAMPLE" =~ ^[A-Za-z0-9_]+$ && "$TARGET" =~ ^[0-9]+$ && "$SEED" =~ ^[0-9]+$ ]] || { echo "Invalid arguments" >&2; exit 2; }
INPUT_DIR="$BASE/variant_calling_pilot_inputs/$SAMPLE"
SAMPLE_TSV="$INPUT_DIR/sample.tsv"
ROUND1_OUT="$BASE/variant_calling_pilot/ds${TARGET}x_seed${SEED}/round1"
ROUND2_OUT="$BASE/variant_calling_pilot/ds${TARGET}x_seed${SEED}/round2"
WORK="$BASE/nf_work/ds${TARGET}x_seed${SEED}/round2"
[[ -s "$SAMPLE_TSV" ]] || { echo "Missing Round1 pilot sample TSV" >&2; exit 1; }
# Guard that the required Round1 data exist; directories alone don't count.
R1_BAM="$ROUND1_OUT/$SAMPLE/round_1/candidate_reads/$SAMPLE.with_mates.bam"
R1_NUMT="$ROUND1_OUT/$SAMPLE/round_1/numt_decoy_ref/$SAMPLE.original_numt.fa"
R1_NUMT_VCF="$ROUND1_OUT/$SAMPLE/round_1/numt_decoy_variant_calling/$SAMPLE.numt_decoy.pass.split.vcf.gz"
R1_VCF_DIR="$ROUND1_OUT/$SAMPLE/round_1/mtdna_variant_calling"
for f in "$R1_BAM" "$R1_NUMT" "$R1_NUMT_VCF"; do
    [[ -s "$f" ]] || { echo "ERROR: Missing Round1 prerequisite $f" >&2; exit 1; }
done
[[ -d "$R1_VCF_DIR" ]] || { echo "ERROR: Missing Round1 mtDNA variant calling: $R1_VCF_DIR" >&2; exit 1; }
REF_MT="${HOM_MT_REF_DIR:-/home/lt692/scratch_pi_njl27/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_chrM}"
CONSENSUS_MIN_DP="${HOM_CONSENSUS_MIN_DP:-100}"
mkdir -p "$ROUND2_OUT" "$WORK" "$BASE/logs"
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "Submitting Round2 HOM pilot sample=$SAMPLE target=$TARGET seed=$SEED"
    sbatch --parsable --job-name="hom_r2_${SAMPLE}" \
      --output="$BASE/logs/hom_r2_${SAMPLE}_${TARGET}x_${SEED}_%j.out" \
      --error="$BASE/logs/hom_r2_${SAMPLE}_${TARGET}x_${SEED}_%j.err" \
      --cpus-per-task=2 --mem=8G --time=24:00:00 \
      --export=ALL,HOM_ANALYSIS_ROOT="$ANALYSIS_ROOT",HOM_MT_REF_DIR="$REF_MT",HOM_CONSENSUS_MIN_DP="$CONSENSUS_MIN_DP" \
      "$0" "$SAMPLE" "$TARGET" "$SEED"
    exit 0
fi
module load Nextflow/24.10.2
cd "$CALL_REPO"
nextflow run "$CALL_REPO/primate_pipeline_round2_consensus_NUMT.nf" \
    -profile cluster -resume -w "$WORK" \
    --sample_tsv "$SAMPLE_TSV" \
    --outdir "$ROUND2_OUT" \
    --round1_outdir "$ROUND1_OUT" \
    --ref_dir "$REF_MT" \
    --consensus_min_dp "$CONSENSUS_MIN_DP"
