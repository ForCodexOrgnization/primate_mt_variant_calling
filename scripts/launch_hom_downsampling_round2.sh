#!/usr/bin/env bash
# Run authoritative Round2 on one completed Round1 HOM downsampling pilot.
# bash scripts/launch_hom_downsampling_round2.sh ERS12091908 100 101
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
SAMPLE="${1:?Provide sample ID}"; TARGET="${2:?Provide raw depth target}"; SEED="${3:?Provide seed}"
case "$TARGET" in 30) DEFAULT_DP=10;; 50) DEFAULT_DP=20;; 100) DEFAULT_DP=30;; *) DEFAULT_DP=100;; esac
CONSENSUS_MIN_DP="${HOM_CONSENSUS_MIN_DP:-$DEFAULT_DP}"
CONSENSUS_MIN_AF="${HOM_CONSENSUS_MIN_AF:-0.95}"
[[ "$CONSENSUS_MIN_DP" =~ ^[1-9][0-9]*$ && "$CONSENSUS_MIN_AF" =~ ^(0\.[0-9]+|1(\.0+)?)$ ]] || { echo "Invalid consensus thresholds" >&2; exit 2; }
CONSENSUS_EXPR="FMT/AF>=${CONSENSUS_MIN_AF} && FMT/DP>=${CONSENSUS_MIN_DP}"
[[ "$SAMPLE" =~ ^[A-Za-z0-9_]+$ && "$TARGET" =~ ^[0-9]+$ && "$SEED" =~ ^[0-9]+$ ]] || { echo "Invalid arguments" >&2; exit 2; }
INPUT_DIR="$BASE/variant_calling_pilot_inputs/$SAMPLE/ds${TARGET}x_seed${SEED}_cdp${CONSENSUS_MIN_DP}"
SAMPLE_TSV="$INPUT_DIR/sample.tsv"
ROUND1_OUT="$BASE/variant_calling_pilot/ds${TARGET}x_seed${SEED}_cdp${CONSENSUS_MIN_DP}/round1"
ROUND2_OUT="$BASE/variant_calling_pilot/ds${TARGET}x_seed${SEED}_cdp${CONSENSUS_MIN_DP}/round2"
WORK="$BASE/nf_work/ds${TARGET}x_seed${SEED}_cdp${CONSENSUS_MIN_DP}/$SAMPLE/round2"
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
REF_MT="${HOM_MT_REF_DIR:-/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_chrM}"
REF_WHOLE="${HOM_GLOBAL_REF_DIR:-/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_whole}"
REF_NAME="$(awk -F '\t' 'NR==1 {print $3; exit}' "$SAMPLE_TSV")"
[[ -s "$REF_WHOLE/$REF_NAME.fa" && -s "$REF_WHOLE/$REF_NAME.fa.fai" ]] || {
    echo "ERROR: Missing matching whole-genome reference/index: $REF_WHOLE/$REF_NAME.fa" >&2
    exit 1
}
CRAM_LINK_ROOT="$INPUT_DIR/cram"
[[ -s "$CRAM_LINK_ROOT/$SAMPLE/alignment/$SAMPLE.cram" && -s "$CRAM_LINK_ROOT/$SAMPLE/alignment/$SAMPLE.cram.crai" ]] || {
    echo "ERROR: Round1 pilot original CRAM links missing: $CRAM_LINK_ROOT/$SAMPLE/alignment" >&2
    exit 1
}
mkdir -p "$ROUND2_OUT" "$WORK" "$BASE/logs"
printf "sample\ttarget_depth\tseed\tconsensus_min_af\tconsensus_min_dp\tconsensus_filter_expr\n%s\t%s\t%s\t%s\t%s\t%s\n" "$SAMPLE" "$TARGET" "$SEED" "$CONSENSUS_MIN_AF" "$CONSENSUS_MIN_DP" "$CONSENSUS_EXPR" > "$ROUND2_OUT/pilot_consensus_parameters.tsv"
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "Submitting Round2 HOM pilot sample=$SAMPLE target=$TARGET seed=$SEED"
    sbatch --parsable --partition="${HOM_SLURM_PARTITION:-ycga}" --job-name="hom_r2_${SAMPLE}" \
      --output="$BASE/logs/hom_r2_${SAMPLE}_${TARGET}x_${SEED}_%j.out" \
      --error="$BASE/logs/hom_r2_${SAMPLE}_${TARGET}x_${SEED}_%j.err" \
      --cpus-per-task=2 --mem=8G --time=24:00:00 \
      --export=ALL,HOM_CALL_REPO="$CALL_REPO",HOM_ANALYSIS_ROOT="$ANALYSIS_ROOT",HOM_MT_REF_DIR="$REF_MT",HOM_GLOBAL_REF_DIR="$REF_WHOLE",HOM_CONSENSUS_MIN_DP="$CONSENSUS_MIN_DP",HOM_CONSENSUS_MIN_AF="$CONSENSUS_MIN_AF" \
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
# Nextflow accepts a comma-separated configuration list. The optional
# experiment-only file overrides the site-wide executor submission limits.
NF_CONFIG="$CALL_REPO/nextflow_mcc.config"
if [[ -n "${HOM_NF_EXTRA_CONFIG:-}" ]]; then
    [[ -s "$HOM_NF_EXTRA_CONFIG" ]] || { echo "ERROR: Missing HOM_NF_EXTRA_CONFIG=$HOM_NF_EXTRA_CONFIG" >&2; exit 1; }
    NF_CONFIG="$NF_CONFIG,$HOM_NF_EXTRA_CONFIG"
fi
echo "INFO: Nextflow config(s): $NF_CONFIG"
# Nextflow creates .nextflow/history.lock relative to the launch directory.
# Run from this pilot's writable, isolated work directory rather than the
# variant-calling source checkout, which might be read-only on McCleary.
mkdir -p "$WORK/.nextflow"
cd "$WORK"
echo "INFO: Nextflow launch directory: $PWD"
nextflow run "$CALL_REPO/primate_pipeline_round2_consensus_NUMT.nf" \
    -c "$NF_CONFIG" -profile cluster -resume -w "$WORK" \
    --sample_tsv "$SAMPLE_TSV" \
    --outdir "$ROUND2_OUT" \
    --round1_outdir "$ROUND1_OUT" \
    --ref_dir "$REF_MT" \
    --global_ref_dir "$REF_WHOLE" \
    --mtcn_nuclear_ref_dir "$REF_WHOLE" \
    --cram_dirs "$CRAM_LINK_ROOT" \
    --wdl_script_round2 "$CALL_REPO/MitochondriaPipeline_multisample_autoLen_round2_NUMT.wdl" \
    --consensus_min_dp "$CONSENSUS_MIN_DP" \
    --consensus_min_af "$CONSENSUS_MIN_AF" \
    --consensus_filter_expr "$CONSENSUS_EXPR"
