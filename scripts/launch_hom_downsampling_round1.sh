#!/usr/bin/env bash
# Round 1 HOM-depth pilot: consume candidate BAMs made by primate_mtDNA_analysis.
# Usage: bash scripts/launch_hom_downsampling_round1.sh SAMPLE DEPTH SEED
#        bash scripts/launch_hom_downsampling_round1.sh ERS12091908 100 101
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# sbatch executes a spool copy of this script. Carry the real checkout path
# from the submit host; BASH_SOURCE on a worker is under /var/spool/slurmd.
CALL_REPO="${HOM_CALL_REPO:-$(cd "$SCRIPT_DIR/.." && pwd)}"
[[ -s "$CALL_REPO/primate_pipeline_numt_decoy_round1.nf" ]] || {
    echo "ERROR: Invalid HOM_CALL_REPO=$CALL_REPO (Round1 workflow missing)" >&2
    exit 1
}
ANALYSIS_ROOT="${HOM_ANALYSIS_ROOT:-/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only}"
BASE="$ANALYSIS_ROOT/results/qc_coverage_only/downsampling_pilot"
SAMPLE="${1:?Provide sample ID}"
TARGET="${2:?Provide raw candidate depth target}"
SEED="${3:?Provide subsampling seed}"
# Explicit experimental condition; optional HOM_CONSENSUS_MIN_DP overrides the adaptive default.
case "$TARGET" in 30) DEFAULT_DP=10;; 50) DEFAULT_DP=20;; 100) DEFAULT_DP=30;; *) DEFAULT_DP=100;; esac
CONSENSUS_MIN_DP="${HOM_CONSENSUS_MIN_DP:-$DEFAULT_DP}"
CONSENSUS_MIN_AF="${HOM_CONSENSUS_MIN_AF:-0.95}"
[[ "$CONSENSUS_MIN_DP" =~ ^[1-9][0-9]*$ && "$CONSENSUS_MIN_AF" =~ ^(0\.[0-9]+|1(\.0+)?)$ ]] || { echo "Invalid consensus thresholds" >&2; exit 2; }
CONSENSUS_EXPR="FMT/AF>=${CONSENSUS_MIN_AF} && FMT/DP>=${CONSENSUS_MIN_DP}"
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
INPUT_DIR="$BASE/variant_calling_pilot_inputs/$SAMPLE/ds${TARGET}x_seed${SEED}_cdp${CONSENSUS_MIN_DP}"
OUT="$BASE/variant_calling_pilot/ds${TARGET}x_seed${SEED}_cdp${CONSENSUS_MIN_DP}/round1"
WORK="$BASE/nf_work/ds${TARGET}x_seed${SEED}_cdp${CONSENSUS_MIN_DP}/$SAMPLE/round1"
CRAM_SRC="${HOM_WGS_ROOT:-/home/lt692/ycga_work/primate_downsampling_pilot/wgs_input}/$SAMPLE/$SAMPLE.cram"
[[ -s "$CRAM_SRC" && -s "$CRAM_SRC.crai" ]] || { echo "ERROR: Original CRAM/CRAI not found at $CRAM_SRC" >&2; exit 1; }
CRAM_LINK_DIR="$INPUT_DIR/cram/$SAMPLE/alignment"
mkdir -p "$CRAM_LINK_DIR" "$OUT" "$WORK" "$BASE/logs"
printf "sample\ttarget_depth\tseed\tconsensus_min_af\tconsensus_min_dp\tconsensus_filter_expr\n%s\t%s\t%s\t%s\t%s\t%s\n" "$SAMPLE" "$TARGET" "$SEED" "$CONSENSUS_MIN_AF" "$CONSENSUS_MIN_DP" "$CONSENSUS_EXPR" > "$OUT/pilot_consensus_parameters.tsv"
# Round2 FIND_ROUND1_OUTPUTS expects the published candidate BAM in Round1
# even when Round1 EXTRACT_CANDIDATE_READS is intentionally bypassed.
ROUND1_CANDIDATE="$OUT/$SAMPLE/round_1/candidate_reads"
mkdir -p "$ROUND1_CANDIDATE"
ln -sfn "$BAM" "$ROUND1_CANDIDATE/$SAMPLE.with_mates.bam"
ln -sfn "$BAM.bai" "$ROUND1_CANDIDATE/$SAMPLE.with_mates.bam.bai"
ln -sfn "$CRAM_SRC" "$CRAM_LINK_DIR/$SAMPLE.cram"
ln -sfn "$CRAM_SRC.crai" "$CRAM_LINK_DIR/$SAMPLE.cram.crai"
SAMPLE_TSV="$INPUT_DIR/sample.tsv"
printf '%s\t%s\t%s\n' "$SAMPLE" "$SPECIES" "$REF_NAME" > "$SAMPLE_TSV"
# Prefer the original, sample-specific production NUMT BED if available.
# For pilots transferred from a different filesystem, reconstruct ONLY the
# nuclear intervals from the exact chrM+NUMT extraction BED. Do not include chrM
# as a NUMT contig: PREPARE_DECOY_REFERENCE adds chrM separately.
NUMT_BED_DIR="${HOM_NUMT_BED_DIR:-/nfs/roberts/pi/pi_njl27/lt692/numt_discovery/primate_results_numt_besthit_strict}"
NUMT_BED="$NUMT_BED_DIR/$SAMPLE.highconf_numt.bed"
REF_WHOLE="${HOM_GLOBAL_REF_DIR:-/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_whole}"
REF_MT="${HOM_MT_REF_DIR:-/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_chrM}"
[[ -s "$REF_WHOLE/$REF_NAME.fa" && -s "$REF_MT/$REF_NAME.fa" ]] || { echo "ERROR: Missing matching .fa reference for $REF_NAME; whole=$REF_WHOLE/$REF_NAME.fa mt=$REF_MT/$REF_NAME.fa" >&2; exit 1; }
if [[ ! -e "$NUMT_BED" ]]; then
    SOURCE_BED="${HOM_ORIGINAL_TARGET_BED:-/home/lt692/ycga_work/primate_downsampling_pilot/candidate_input/$SAMPLE/round_1/candidate_reads/$SAMPLE.target_chrM_plus_NUMT.bed}"
    [[ -s "$SOURCE_BED" ]] || { echo "ERROR: Neither production NUMT BED nor extraction BED found: $NUMT_BED / $SOURCE_BED" >&2; exit 1; }
    NUMT_BED_DIR="$INPUT_DIR/numt_bed_from_extraction"
    mkdir -p "$NUMT_BED_DIR"
    NUMT_BED="$NUMT_BED_DIR/$SAMPLE.highconf_numt.bed"
    TMP_BED="$NUMT_BED.tmp.$"
    awk 'BEGIN{FS=OFS="\t"}
      /^[#]/ || NF==0 {next}
      NF<3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3<=$2 {bad++;next}
      $1=="chrM" {next}
      {print $1,$2,$3; kept++}
      END{if(bad>0 || kept==0){print "ERROR: Invalid or missing nuclear NUMT intervals in extraction BED" > "/dev/stderr";exit 1}}
    ' "$SOURCE_BED" > "$TMP_BED" || { rm -f "$TMP_BED"; exit 1; }
    # Fail fast if the extraction interval set does not match the selected
    # whole-genome FASTA coordinate system.
    [[ -s "$REF_WHOLE/$REF_NAME.fa.fai" ]] || { echo "ERROR: Missing whole-genome FASTA index" >&2; rm -f "$TMP_BED"; exit 1; }
    awk 'BEGIN{FS=OFS="\t";bad=0}
      NR==FNR {lens[$1]=$2;next}
      !($1 in lens) || $3>lens[$1] {print "ERROR: NUMT BED interval not present in reference: " $0 > "/dev/stderr";bad++}
      END{exit (bad>0)}
    ' "$REF_WHOLE/$REF_NAME.fa.fai" "$TMP_BED" || { rm -f "$TMP_BED"; exit 1; }
    mv "$TMP_BED" "$NUMT_BED"
    echo "INFO: Production NUMT BED unavailable; derived $(wc -l < "$NUMT_BED") nuclear intervals from original extraction BED: $SOURCE_BED"
    echo "WARNING: Verify the derived intervals match the original Round1 decoy BED before interpreting HOM results."
fi
echo "INFO: Round1 NUMT BED: $NUMT_BED"
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "Submitting Round1 HOM pilot sample=$SAMPLE target=$TARGET seed=$SEED; consensus_min_dp=$CONSENSUS_MIN_DP"
    sbatch --parsable --partition="${HOM_SLURM_PARTITION:-ycga}" --job-name="hom_r1_${SAMPLE}" \
      --output="$BASE/logs/hom_r1_${SAMPLE}_${TARGET}x_${SEED}_%j.out" \
      --error="$BASE/logs/hom_r1_${SAMPLE}_${TARGET}x_${SEED}_%j.err" \
      --cpus-per-task=2 --mem=8G --time=24:00:00 \
      --export=ALL,HOM_CALL_REPO="$CALL_REPO",HOM_ANALYSIS_ROOT="$ANALYSIS_ROOT",HOM_REF_NAME="$REF_NAME",HOM_CONSENSUS_MIN_DP="$CONSENSUS_MIN_DP",HOM_CONSENSUS_MIN_AF="$CONSENSUS_MIN_AF",HOM_NUMT_BED_DIR="$NUMT_BED_DIR",HOM_GLOBAL_REF_DIR="$REF_WHOLE",HOM_MT_REF_DIR="$REF_MT" \
      "$0" "$SAMPLE" "$TARGET" "$SEED"
    exit 0
fi
load_pilot_nextflow() {
    # Prefer an existing executable (including user conda installations).
    if command -v nextflow >/dev/null 2>&1; then
        echo "INFO: Using Nextflow: $(command -v nextflow)"
        nextflow -version
        return
    fi
    if [[ -n "${NEXTFLOW_MODULE:-}" ]]; then
        module load "$NEXTFLOW_MODULE"
        command -v nextflow >/dev/null 2>&1 || {
            echo "ERROR: NEXTFLOW_MODULE=$NEXTFLOW_MODULE loaded but nextflow is not on PATH" >&2
            exit 1
        }
        echo "INFO: Using Nextflow module $NEXTFLOW_MODULE"
        nextflow -version
        return
    fi
    echo "ERROR: Nextflow not found on PATH. On McCleary, check: module spider Nextflow" >&2
    echo "Then either activate an environment providing nextflow or set NEXTFLOW_MODULE to a valid module before submitting." >&2
    exit 1
}
load_pilot_nextflow
# Nextflow creates .nextflow/history.lock relative to the launch directory.
# Run from this pilot's writable, isolated work directory rather than the
# variant-calling source checkout, which might be read-only on McCleary.
mkdir -p "$WORK/.nextflow"
cd "$WORK"
echo "INFO: Nextflow launch directory: $PWD"
nextflow run "$CALL_REPO/primate_pipeline_numt_decoy_round1.nf" \
    -c "$CALL_REPO/nextflow_mcc.config" -profile cluster -resume -w "$WORK" \
    --sample_tsv "$SAMPLE_TSV" \
    --outdir "$OUT" \
    --cram_dirs "$INPUT_DIR/cram" \
    --numt_bed_dir "$NUMT_BED_DIR" \
    --global_ref_dir "$REF_WHOLE" \
    --ref_dir "$REF_MT" \
    --wdl_script "$CALL_REPO/MitochondriaPipeline_multisample_autoLen.wdl" \
    --candidate_bam_root "$BAM_ROOT" \
    --candidate_target_depth "$TARGET" \
    --candidate_seed "$SEED" \
    --consensus_min_dp "$CONSENSUS_MIN_DP" \
    --consensus_min_af "$CONSENSUS_MIN_AF" \
    --consensus_filter_expr "$CONSENSUS_EXPR"
