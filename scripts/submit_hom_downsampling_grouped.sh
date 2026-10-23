#!/usr/bin/env bash
# Multi-sample HOM pilot recovery. One Slurm array element = 10 samples at
# identical depth/seed/consensus-DP. Each element runs one Round1 Nextflow,
# then one Round2 Nextflow, preserving existing published output directories.
#
# Examples:
#   HOM_CONSENSUS_MIN_DP=100 bash scripts/submit_hom_downsampling_grouped.sh --dry-run
#   HOM_CONSENSUS_MIN_DP=100 bash scripts/submit_hom_downsampling_grouped.sh --array 1
#   unset HOM_CONSENSUS_MIN_DP
#   bash scripts/submit_hom_downsampling_grouped.sh
set -euo pipefail
SELF="${BASH_SOURCE[0]}"
SCRIPT_DIR="$(cd "$(dirname "$SELF")" && pwd)"
CALL_REPO="${HOM_CALL_REPO:-$(cd "$SCRIPT_DIR/.." && pwd)}"
ANALYSIS_ROOT="${HOM_ANALYSIS_ROOT:-/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only}"
BASE="$ANALYSIS_ROOT/results/qc_coverage_only/downsampling_pilot"
PARTITION="${HOM_SLURM_PARTITION:-ycga}"
NEXTFLOW_MODULE="${NEXTFLOW_MODULE:-Nextflow/24.04.2}"
RATE="${HOM_GROUP_SUBMIT_RATE:-120/h}"
QUEUE="${HOM_GROUP_QUEUE_SIZE:-10}"
GAP="${HOM_GROUP_SUCCESS_GAP_SECONDS:-1200}"
[[ "$RATE" =~ ^[1-9][0-9]*/h$ && "$QUEUE" =~ ^[1-9][0-9]*$ && "$GAP" =~ ^[0-9]+$ ]] || {
 echo "Invalid RATE/QUEUE/GAP" >&2; exit 2;
}
REF_WHOLE="${HOM_GLOBAL_REF_DIR:-/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_whole}"
REF_MT="${HOM_MT_REF_DIR:-/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_chrM}"
WGS_ROOT="${HOM_WGS_ROOT:-/home/lt692/ycga_work/primate_downsampling_pilot/wgs_input}"
META_ROOT="${HOM_CANDIDATE_META_ROOT:-/home/lt692/ycga_work/primate_downsampling_pilot/candidate_input}"
PROD_BED="${HOM_NUMT_BED_DIR:-/nfs/roberts/pi/pi_njl27/lt692/numt_discovery/primate_results_numt_besthit_strict}"
CONS_AF="${HOM_CONSENSUS_MIN_AF:-0.95}"

load_nextflow() {
 if ! command -v nextflow >/dev/null 2>&1; then
   module load "$NEXTFLOW_MODULE"
 fi
 command -v nextflow >/dev/null 2>&1 || { echo "Nextflow not available" >&2; exit 1; }
}
if [[ "${1:-}" == "--worker" ]]; then
 INDEX="${SLURM_ARRAY_TASK_ID:?Must execute as Slurm array task}"
 case "$INDEX" in
   1) DEPTH=30; SEED=101 ;;
   2) DEPTH=30; SEED=202 ;;
   3) DEPTH=50; SEED=101 ;;
   4) DEPTH=50; SEED=202 ;;
   5) DEPTH=100; SEED=101 ;;
   6) DEPTH=100; SEED=202 ;;
   *) echo "Invalid group index: $INDEX" >&2; exit 2 ;;
 esac
 case "$DEPTH" in 30) ADAPT_DP=10;; 50) ADAPT_DP=20;; 100) ADAPT_DP=30;; esac
 CDP="${HOM_CONSENSUS_MIN_DP:-$ADAPT_DP}"
 TAG="ds${DEPTH}x_seed${SEED}_cdp${CDP}"
 GROUP="$BASE/variant_calling_pilot_group_inputs/$TAG"
 OUT="$BASE/variant_calling_pilot/$TAG"
 WORK="$BASE/nf_work_grouped/$TAG"
 echo "GROUP START=$TAG date=$(date -Is) NF_RATE=$RATE QUEUE=$QUEUE"
 python3 "$CALL_REPO/scripts/run_hom_downsampling_group.py" \
   --base "$BASE" --depth "$DEPTH" --seed "$SEED" --cdp "$CDP" \
   --wgs-root "$WGS_ROOT" --candidate-meta-root "$META_ROOT" \
   --prod-numt-dir "$PROD_BED" --ref-whole "$REF_WHOLE" --ref-mt "$REF_MT"
 mkdir -p "$OUT/round1" "$OUT/round2" "$WORK" "$BASE/logs"
 CONFIG="$GROUP/nextflow_throttle.config"
 # This configuration is scoped to one multi-sample Nextflow invocation.
 cat > "$CONFIG" <<EOF
executor {
    queueSize = $QUEUE
    submitRateLimit = '$RATE'
}
EOF
 NF_CONFIG="$CALL_REPO/nextflow_mcc.config,$CONFIG"
 SAMPLES="$GROUP/sample.tsv"
 R1_PENDING=0
 R2_PENDING=0
 while IFS=$'\t' read -r SAMPLE SPECIES REF; do
   [[ -n "$SAMPLE" ]] || continue
   if ! bash "$CALL_REPO/scripts/round1_outputs_complete.sh" "$OUT/round1" "$SAMPLE" >/dev/null 2>&1; then
     R1_PENDING=$(( R1_PENDING + 1 ))
   fi
   ROOT2="$OUT/round2/$SAMPLE"
   VCF="$ROOT2/round_2_variant_calling_original_coords/$SAMPLE.round2.original_coords.clean.final.split.vcf.gz"
   COV="$ROOT2/round_2_variant_calling_original_coords/$SAMPLE.round2.original_coords.per_base_coverage.tsv"
   MTCN="$ROOT2/round_2/mtcn/$SAMPLE.round2.mtcn.tsv"
   [[ -s "$VCF" && -s "$VCF.tbi" && -s "$COV" && -s "$MTCN" ]] || R2_PENDING=$(( R2_PENDING + 1 ))
 done < "$SAMPLES"
 echo "STATUS: round1_pending=$R1_PENDING round2_pending=$R2_PENDING"
 if (( R1_PENDING == 0 && R2_PENDING == 0 )); then
   echo "SKIP GROUP all ten samples completed"
   exit 0
 fi
 load_nextflow
 EXPR="FMT/AF>=${CONS_AF} && FMT/DP>=${CDP}"
 # Round1 uses the original production multi-sample workflow, including its
 # per-sample completed-output detection.
 if (( R1_PENDING > 0 )); then
   mkdir -p "$WORK/round1/.nextflow"
   cd "$WORK/round1"
   nextflow run "$CALL_REPO/primate_pipeline_numt_decoy_round1.nf" \
     -c "$NF_CONFIG" -profile cluster -resume -w "$WORK/round1" \
     --sample_tsv "$SAMPLES" --outdir "$OUT/round1" \
     --cram_dirs "$GROUP/cram" --numt_bed_dir "$GROUP/numt_bed" \
     --global_ref_dir "$REF_WHOLE" --ref_dir "$REF_MT" \
     --wdl_script "$CALL_REPO/MitochondriaPipeline_multisample_autoLen.wdl" \
     --candidate_bam_root "$BASE/candidate_downsampled" \
     --candidate_target_depth "$DEPTH" --candidate_seed "$SEED" \
     --consensus_min_dp "$CDP" --consensus_min_af "$CONS_AF" \
     --consensus_filter_expr "$EXPR"
 fi
 # Verify that the Round1 outputs satisfy the production contract before
 # giving all ten samples to the next workflow.
 while IFS=$'\t' read -r SAMPLE SPECIES REF; do
   [[ -n "$SAMPLE" ]] || continue
   bash "$CALL_REPO/scripts/round1_outputs_complete.sh" "$OUT/round1" "$SAMPLE" || {
      echo "ERROR: Round1 incomplete after group execution: $SAMPLE" >&2; exit 1;
   }
 done < "$SAMPLES"
 if (( R2_PENDING > 0 )); then
   mkdir -p "$WORK/round2/.nextflow"
   cd "$WORK/round2"
   nextflow run "$CALL_REPO/primate_pipeline_round2_consensus_NUMT.nf" \
     -c "$NF_CONFIG" -profile cluster -resume -w "$WORK/round2" \
     --sample_tsv "$SAMPLES" --outdir "$OUT/round2" \
     --round1_outdir "$OUT/round1" --ref_dir "$REF_MT" \
     --global_ref_dir "$REF_WHOLE" --mtcn_nuclear_ref_dir "$REF_WHOLE" \
     --cram_dirs "$GROUP/cram" \
     --wdl_script_round2 "$CALL_REPO/MitochondriaPipeline_multisample_autoLen_round2_NUMT.wdl" \
     --consensus_min_dp "$CDP" --consensus_min_af "$CONS_AF" \
     --consensus_filter_expr "$EXPR"
 fi
 while IFS=$'\t' read -r SAMPLE SPECIES REF; do
   [[ -n "$SAMPLE" ]] || continue
   R="$OUT/round2/$SAMPLE"
   F="$R/round_2_variant_calling_original_coords/$SAMPLE.round2.original_coords"
   [[ -s "$F.clean.final.split.vcf.gz" && -s "$F.clean.final.split.vcf.gz.tbi" && -s "$F.per_base_coverage.tsv" && -s "$R/round_2/mtcn/$SAMPLE.round2.mtcn.tsv" ]] || {
     echo "ERROR: Round2 incomplete after group execution: $SAMPLE" >&2; exit 1;
   }
 done < "$SAMPLES"
 echo "SUCCESS GROUP $TAG date=$(date -Is)"
 # Delay only between groups that actually ran Nextflow (never skipped groups).
 # This reduces, but cannot enforce, the account-wide rolling-hour quota.
 if (( GAP > 0 )); then
   echo "GROUP_GAP ${GAP}s"
   sleep "$GAP"
 fi
 exit 0
fi
DRY=0
RANGE="1-6"
while (( $# > 0 )); do
 case "$1" in
  --dry-run) DRY=1; shift ;;
  --array) [[ $# -ge 2 ]] || exit 2; RANGE="$2"; shift 2 ;;
  *) echo "Usage: $0 [--dry-run] [--array 1-6|INDEX]" >&2; exit 2 ;;
 esac
done
[[ "$RANGE" =~ ^[1-6](-[1-6])?$ ]] || { echo "Invalid array range $RANGE" >&2; exit 2; }
mkdir -p "$BASE/logs"
echo "Grouped HOM recovery: array=$RANGE max_groups=1 sample_count_per_group=10"
echo "Consensus DP=${HOM_CONSENSUS_MIN_DP:-adaptive} AF=$CONS_AF; Nextflow rate=$RATE queue=$QUEUE gap=$GAP"
echo "IMPORTANT: Do not overlap this job with the older per-sample recovery array."
if (( DRY )); then
 echo "DRY RUN: no jobs submitted."
 exit 0
fi
# Job array is a single scheduler submission; Nextflow submits its internal
# processes subject to RATE within the group.
JOB="$(sbatch --parsable --partition="$PARTITION" --array="${RANGE}%1" \
  --job-name=hom_grouped --cpus-per-task=2 --mem=8G --time=24:00:00 \
  --output="$BASE/logs/hom_grouped_%A_%a.out" \
  --error="$BASE/logs/hom_grouped_%A_%a.err" \
  --export=ALL,HOM_CALL_REPO="$CALL_REPO",HOM_ANALYSIS_ROOT="$ANALYSIS_ROOT",NEXTFLOW_MODULE="$NEXTFLOW_MODULE" \
  "$CALL_REPO/scripts/submit_hom_downsampling_grouped.sh" --worker)"
JOB="${JOB%%;*}"
printf 'job_id\tmode\tarray\trate\tqueue\n%s\t%s\t%s\t%s\t%s\n' \
 "$JOB" "${HOM_CONSENSUS_MIN_DP:-adaptive}" "$RANGE" "$RATE" "$QUEUE" >> "$BASE/hom_grouped_jobs.tsv"
echo "Submitted grouped array $JOB"
