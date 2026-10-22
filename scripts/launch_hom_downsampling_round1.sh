#!/usr/bin/env bash
# Round 1 HOM-depth pilot: consume candidate BAMs made by primate_mtDNA_analysis.
# Usage: bash scripts/launch_hom_downsampling_round1.sh SAMPLE DEPTH SEED
#        bash scripts/launch_hom_downsampling_round1.sh ERS12091908 100 101
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CALL_REPO="$(cd "$SCRIPT_DIR/.." && pwd)"
ANALYSIS_ROOT="${HOM_ANALYSIS_ROOT:-/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only}"
BASE="$ANALYSIS_ROOT/results/qc_coverage_only/downsampling_pilot"
SAMPLE="${1:?Provide sample ID}"
TARGET="${2:?Provide raw candidate depth target}"
SEED="${3:?Provide subsampling seed}"
[[ "$SAMPLE" =~ ^[A-Za-z0-9_]+$ && "$TARGET" =~ ^[0-9]+$ && "$SEED" =~ ^[0-9]+$ ]] || { echo "Invalid sample/target/seed" >&2; exit 2; }
BAM_ROOT="$BASE/candidate_downsampled"
BAM="$BAM_ROOT/$SAMPLE/$SAMPLE.candidate.ds${TARGET}x.seed${SEED}.bam"
[[ -s "$BAM" && -s "$BAM.bai" ]] || { echo "ERROR: Missing BAM/index: $BAM" >&2; exit 1; }
PILOT_TABLE="$BASE/pilot_samples.tsv"
[[ -s "$PILOT_TABLE" ]] || { echo "Missing sample table: $PILOT_TABLE" >&2; exit 1; }
SPECIES="$(awk -F '\t' -v s="$SAMPLE" 'NR>1 && $1==s {print $2;exit}' "$PILOT_TABLE")"
[[ -n "$SPECIES" ]] || { echo "Sample $SAMPLE not present in $PILOT_TABLE" >&2; exit 1; }
# The pilot list uses species-name references. Override for any exceptions.
REF_NAME="${HOM_REF_NAME:-$SPECIES}"
INPUT_DIR="$BASE/variant_calling_pilot_inputs/$SAMPLE"
OUT="$BASE/variant_calling_pilot/ds${TARGET}x_seed${SEED}/round1"
WORK="$BASE/nf_work/ds${TARGET}x_seed${SEED}/round1"
CRAM_SRC="${HOM_WGS_ROOT:-/home/lt692/ycga_work/primate_downsampling_pilot/wgs_input}/$SAMPLE/$SAMPLE.cram"
[[ -s "$CRAM_SRC" && -s "$CRAM_SRC.crai" ]] || { echo "ERROR: Original CRAM/CRAI not found at $CRAM_SRC" >&2; exit 1; }
CRAM_LINK_DIR="$INPUT_DIR/cram/$SAMPLE/alignment"
mkdir -p "$CRAM_LINK_DIR" "$OUT" "$WORK" "$BASE/logs"
ln -sfn "$CRAM_SRC" "$CRAM_LINK_DIR/$SAMPLE.cram"
ln -sfn "$CRAM_SRC.crai" "$CRAM_LINK_DIR/$SAMPLE.cram.crai"
SAMPLE_TSV="$INPUT_DIR/sample.tsv"
printf '%s\t%s\t%s\n' "$SAMPLE" "$SPECIES" "$REF_NAME" > "$SAMPLE_TSV"
# Reuse production NUMT BEDs and reference config by default. Do not derive
# NUMT intervals from the read-extraction BED, which can have different semantics.
NUMT_BED_DIR="${HOM_NUMT_BED_DIR:-/nfs/roberts/pi/pi_njl27/lt692/numt_discovery/primate_results_numt_besthit_strict}"
NUMT_BED="$NUMT_BED_DIR/$SAMPLE.highconf_numt.bed"
[[ -e "$NUMT_BED" ]] || { echo "ERROR: NUMT BED missing: $NUMT_BED" >&2; exit 1; }
REF_WHOLE="${HOM_GLOBAL_REF_DIR:-/home/lt692/scratch_pi_njl27/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_whole}"
REF_MT="${HOM_MT_REF_DIR:-/home/lt692/scratch_pi_njl27/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_chrM}"
[[ -s "$REF_WHOLE/$REF_NAME.fa" && -s "$REF_MT/$REF_NAME.fa" ]] || { echo "Missing matching .fa reference for $REF_NAME" >&2; exit 1; }
# The consensus threshold is not silently changed for the low-depth study.
CONSENSUS_MIN_DP="${HOM_CONSENSUS_MIN_DP:-100}"
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "Submitting Round1 HOM pilot sample=$SAMPLE target=$TARGET seed=$SEED; consensus_min_dp=$CONSENSUS_MIN_DP"
    sbatch --parsable --job-name="hom_r1_${SAMPLE}" \
      --output="$BASE/logs/hom_r1_${SAMPLE}_${TARGET}x_${SEED}_%j.out" \
      --error="$BASE/logs/hom_r1_${SAMPLE}_${TARGET}x_${SEED}_%j.err" \
      --cpus-per-task=2 --mem=8G --time=24:00:00 \
      --export=ALL,HOM_ANALYSIS_ROOT="$ANALYSIS_ROOT",HOM_REF_NAME="$REF_NAME",HOM_CONSENSUS_MIN_DP="$CONSENSUS_MIN_DP",HOM_NUMT_BED_DIR="$NUMT_BED_DIR",HOM_GLOBAL_REF_DIR="$REF_WHOLE",HOM_MT_REF_DIR="$REF_MT" \
      "$0" "$SAMPLE" "$TARGET" "$SEED"
    exit 0
fi
module load Nextflow/24.10.2
cd "$CALL_REPO"
nextflow run "$CALL_REPO/primate_pipeline_numt_decoy_round1.nf" \
    -profile cluster -resume -w "$WORK" \
    --sample_tsv "$SAMPLE_TSV" \
    --outdir "$OUT" \
    --cram_dirs "$INPUT_DIR/cram" \
    --numt_bed_dir "$NUMT_BED_DIR" \
    --global_ref_dir "$REF_WHOLE" \
    --ref_dir "$REF_MT" \
    --candidate_bam_root "$BAM_ROOT" \
    --candidate_target_depth "$TARGET" \
    --candidate_seed "$SEED" \
    --consensus_min_dp "$CONSENSUS_MIN_DP"
