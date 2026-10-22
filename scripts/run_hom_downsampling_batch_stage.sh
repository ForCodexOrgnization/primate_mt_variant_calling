#!/usr/bin/env bash
# Invoked as a Slurm array task by submit_hom_downsampling_batch.sh.
# stage: round1 | round2. Each array index matches the candidate manifest row.
set -euo pipefail
STAGE="${1:?Specify round1 or round2}"
[[ "$STAGE" == "round1" || "$STAGE" == "round2" ]] || { echo "Invalid stage $STAGE" >&2; exit 2; }
TASK="${SLURM_ARRAY_TASK_ID:?Must run as a Slurm array task}"
CALL_REPO="${HOM_CALL_REPO:?Set HOM_CALL_REPO to the full calling checkout}"
ROOT="${HOM_ANALYSIS_ROOT:-/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only}"
BASE="$ROOT/results/qc_coverage_only/downsampling_pilot"
MANIFEST="$BASE/candidate_downsampling_jobs.tsv"
[[ -s "$MANIFEST" ]] || { echo "Missing $MANIFEST" >&2; exit 1; }
# manifest header: sample species input_bam mt_contig original_raw_median target_depth seed fraction
read -r SAMPLE TARGET SEED < <(awk -F '\t' -v n="$TASK" 'NR==n+1 {print $1, $6, $7; exit}' "$MANIFEST")
[[ -n "${SAMPLE:-}" && "$TARGET" =~ ^[0-9]+$ && "$SEED" =~ ^[0-9]+$ ]] || { echo "Invalid manifest task $TASK" >&2; exit 1; }
case "$TARGET" in 30) DEFAULT_DP=10;; 50) DEFAULT_DP=20;; 100) DEFAULT_DP=30;; *) DEFAULT_DP=100;; esac
CONSENSUS_MIN_DP="${HOM_CONSENSUS_MIN_DP:-$DEFAULT_DP}"
echo "HOM_STAGE=$STAGE task=$TASK sample=$SAMPLE target=$TARGET seed=$SEED host=$(hostname)"
R1="$BASE/variant_calling_pilot/ds${TARGET}x_seed${SEED}_cdp${CONSENSUS_MIN_DP}/round1/$SAMPLE"
R2="$BASE/variant_calling_pilot/ds${TARGET}x_seed${SEED}_cdp${CONSENSUS_MIN_DP}/round2/$SAMPLE"
if [[ "$STAGE" == "round1" ]]; then
  # Only skip if both mtDNA split calls and NUMT split calls are present.
  MTDNA="$R1/round_1/mtdna_variant_calling"
  NUMT="$R1/round_1/numt_decoy_variant_calling/$SAMPLE.numt_decoy.pass.split.vcf.gz"
  if [[ -s "$NUMT" && -d "$MTDNA" ]]; then
    SPLIT="$(find "$MTDNA" -type f -name "$SAMPLE.numt_decoy.clean.final.split.vcf" -print -quit)"
    if [[ -n "$SPLIT" && -s "$SPLIT" ]]; then
      echo "SKIP existing completed Round1: $SAMPLE $TARGET $SEED"
      exit 0
    fi
  fi
else
  VCF="$R2/round_2_variant_calling_original_coords/$SAMPLE.round2.original_coords.clean.final.split.vcf.gz"
  COV="$R2/round_2_variant_calling_original_coords/$SAMPLE.round2.original_coords.per_base_coverage.tsv"
  MT="$R2/round_2/mtcn/$SAMPLE.round2.mtcn.tsv"
  if [[ -s "$VCF" && -s "$VCF.tbi" && -s "$COV" && -s "$MT" ]]; then
    echo "SKIP existing completed Round2: $SAMPLE $TARGET $SEED"
    exit 0
  fi
fi
export HOM_CALL_REPO="$CALL_REPO"
export NEXTFLOW_MODULE="${NEXTFLOW_MODULE:-Nextflow/24.04.2}"
bash "$CALL_REPO/scripts/launch_hom_downsampling_${STAGE}.sh" "$SAMPLE" "$TARGET" "$SEED"
