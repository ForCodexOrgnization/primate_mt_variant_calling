#!/usr/bin/env bash
#SBATCH --job-name=NF_Primate_Stream
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=log_streaming/nf_streaming_%A_%a.log

set -euo pipefail

# ==============================================================================
# Per-sample streaming launcher for preprocessing -> NUMT discovery -> round 1 -> round 2.
#
# Run from the repository directory on a login/submission node:
#   bash launch_pipeline_streaming_per_sample.sh
# or submit the coordinator itself:
#   sbatch launch_pipeline_streaming_per_sample.sh
#
# In coordinator mode, this script submits itself as a 1-based Slurm array using:
#   --array=1-N%MAX_CONCURRENT
# Each array task reads one line from FULL_SAMPLE_LIST with sed -n "${SLURM_ARRAY_TASK_ID}p"
# and processes that sample from start to finish without waiting for other samples.
#
# Set CLEAN_ON_SUCCESS=1 to delete only this sample's Nextflow work directories after
# all stages and validations pass. Work directories are always retained on failure.
# ===============================================================================

FULL_SAMPLE_LIST="${FULL_SAMPLE_LIST:-/nfs/roberts/project/pi_njl27/lt692/primate_mt_variant_calling/human_sample.txt}"
PRE_OUTPUT_DIR="${PRE_OUTPUT_DIR:-/nfs/roberts/scratch/pi_njl27/lt692/primate_results_test}"
ROUND_OUTPUT_DIR="${ROUND_OUTPUT_DIR:-/nfs/roberts/scratch/pi_njl27/lt692/primate_results_test}"
ROUND1_OUTDIR="${ROUND1_OUTDIR:-${ROUND_OUTPUT_DIR}}"

NUMT_DISCOVERY_OUTROOT="${NUMT_DISCOVERY_OUTROOT:-${ROUND_OUTPUT_DIR}/numt_discovery}"
NUMT_BESTHIT_OUTDIR="${NUMT_BESTHIT_OUTDIR:-${ROUND_OUTPUT_DIR}/numt_besthit}"
GLOBAL_REF_DIR="${GLOBAL_REF_DIR:-/home/lt692/scratch_pi_njl27/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_whole}"
REF_DIR="${REF_DIR:-/home/lt692/scratch_pi_njl27/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_chrM}"
NUCLEAR_ONLY_REF_DIR="${NUCLEAR_ONLY_REF_DIR:-/home/lt692/scratch_pi_njl27/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_nuclear_only}"

NF_BASE_WORK_DIR="${NF_BASE_WORK_DIR:-/nfs/roberts/scratch/pi_njl27/lt692/nf_work_dir_streaming_per_sample}"
MAX_CONCURRENT="${MAX_CONCURRENT:-10}"
CLEAN_ON_SUCCESS="${CLEAN_ON_SUCCESS:-0}"
LOG_DIR="${LOG_DIR:-log_streaming}"

log() {
    printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" >&2
}

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

validate_file_nonempty() {
    local path="$1"
    [[ -s "${path}" ]] || fail "Missing or empty required file: ${path}"
}

validate_file_exists() {
    local path="$1"
    [[ -e "${path}" ]] || fail "Missing required file: ${path}"
}

# When submitted with sbatch, Slurm may execute a spool copy. Prefer the original
# submission directory when it looks like this repository checkout.
if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/preprocessing.nf" ]]; then
    REPO_DIR="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
else
    REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
cd "${REPO_DIR}"

[[ -f "${FULL_SAMPLE_LIST}" ]] || fail "FULL_SAMPLE_LIST does not exist: ${FULL_SAMPLE_LIST}"
mkdir -p "${LOG_DIR}" "${NF_BASE_WORK_DIR}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    NUM_SAMPLES="$(wc -l < "${FULL_SAMPLE_LIST}" | tr -d '[:space:]')"
    [[ "${NUM_SAMPLES}" -gt 0 ]] || fail "No sample lines found in ${FULL_SAMPLE_LIST}"

    log "Submitting per-sample streaming array for ${NUM_SAMPLES} samples with concurrency ${MAX_CONCURRENT}"
    sbatch --export=ALL --array=1-"${NUM_SAMPLES}"%"${MAX_CONCURRENT}" "$0"
    log "Job array submitted. Monitor with: squeue -u \$USER"
    exit 0
fi

module load Nextflow/24.10.2

log "Running per-sample streaming worker for 1-based task ${SLURM_ARRAY_TASK_ID}"
SAMPLE_LINE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FULL_SAMPLE_LIST}")"
[[ -n "${SAMPLE_LINE}" ]] || fail "No sample line found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

SAMPLE_ID="$(printf '%s\n' "${SAMPLE_LINE}" | awk -F '\t' '{print $1}')"
[[ -n "${SAMPLE_ID}" ]] || fail "Could not extract SAMPLE_ID from line ${SLURM_ARRAY_TASK_ID}"
[[ "${SAMPLE_ID}" != "sample" && "${SAMPLE_ID}" != "sample_id" ]] || fail "Task ${SLURM_ARRAY_TASK_ID} selected a header line (${SAMPLE_ID}); remove headers from FULL_SAMPLE_LIST for 1-based streaming mode"

SAMPLE_WORK_ROOT="${NF_BASE_WORK_DIR}/${SAMPLE_ID}"
SAMPLE_TSV_DIR="${SAMPLE_WORK_ROOT}/inputs"
SAMPLE_TSV="${SAMPLE_TSV_DIR}/${SAMPLE_ID}.sample.tsv"
PRE_WORK_DIR="${SAMPLE_WORK_ROOT}/pre"
ROUND1_WORK_DIR="${SAMPLE_WORK_ROOT}/round1"
ROUND2_WORK_DIR="${SAMPLE_WORK_ROOT}/round2"
# Run each Nextflow invocation from a sample/stage-specific launch directory so
# parallel array tasks do not contend for the repository-level .nextflow cache lock.
PRE_LAUNCH_DIR="${SAMPLE_WORK_ROOT}/nextflow_launch/pre"
ROUND1_LAUNCH_DIR="${SAMPLE_WORK_ROOT}/nextflow_launch/round1"
ROUND2_LAUNCH_DIR="${SAMPLE_WORK_ROOT}/nextflow_launch/round2"
NUMT_CONFIG="${SAMPLE_WORK_ROOT}/${SAMPLE_ID}.numt.config"

mkdir -p "${SAMPLE_TSV_DIR}" "${PRE_WORK_DIR}" "${ROUND1_WORK_DIR}" "${ROUND2_WORK_DIR}" \
    "${PRE_LAUNCH_DIR}" "${ROUND1_LAUNCH_DIR}" "${ROUND2_LAUNCH_DIR}" \
    "${NUMT_DISCOVERY_OUTROOT}" "${NUMT_BESTHIT_OUTDIR}"
printf '%s\n' "${SAMPLE_LINE}" > "${SAMPLE_TSV}"

log "Sample: ${SAMPLE_ID}"
log "One-sample TSV: ${SAMPLE_TSV}"

CRAM_PATH="${PRE_OUTPUT_DIR}/${SAMPLE_ID}/alignment/${SAMPLE_ID}.cram"
CRAI_PATH="${PRE_OUTPUT_DIR}/${SAMPLE_ID}/alignment/${SAMPLE_ID}.cram.crai"

if [[ -s "${CRAM_PATH}" && -s "${CRAI_PATH}" ]]; then
    log "Skipping preprocessing for ${SAMPLE_ID}; existing CRAM and CRAI were found"
else
    log "Starting preprocessing for ${SAMPLE_ID}"
    (
        cd "${PRE_LAUNCH_DIR}"
        nextflow run "${REPO_DIR}/preprocessing.nf" \
            -profile cluster \
            -resume \
            -w "${PRE_WORK_DIR}" \
            --sample_tsv "${SAMPLE_TSV}" \
            --outdir "${PRE_OUTPUT_DIR}"
    )
fi

validate_file_nonempty "${CRAM_PATH}"
validate_file_nonempty "${CRAI_PATH}"
log "Preprocessing validation passed for ${SAMPLE_ID}"

NUMT_BED="${NUMT_BESTHIT_OUTDIR}/${SAMPLE_ID}.highconf_numt.bed"

if [[ -e "${NUMT_BED}" ]]; then
    # Empty high-confidence NUMT BED files are valid for samples with no calls.
    log "Skipping NUMT discovery for ${SAMPLE_ID}; existing high-confidence NUMT BED was found: ${NUMT_BED}"
else
    NUMT_DIR="${REPO_DIR}/numt_detection"
    [[ -d "${NUMT_DIR}" && -f "${NUMT_DIR}/run_numt_end2end.sh" ]] || fail "Bundled numt_detection is missing run_numt_end2end.sh: ${NUMT_DIR}"
    cat > "${NUMT_CONFIG}" <<CONFIG
# Auto-generated by primate_mt_variant_calling/launch_pipeline_streaming_per_sample.sh for ${SAMPLE_ID}
SAMPLES_TSV=${SAMPLE_TSV}
CRAM_ROOT_1=${PRE_OUTPUT_DIR}
CRAM_ROOT_2=${PRE_OUTPUT_DIR}
WHOLE_REF_DIR=${GLOBAL_REF_DIR}
NUCLEAR_ONLY_REF_DIR=${NUCLEAR_ONLY_REF_DIR}
CHRM_REF_DIR=${REF_DIR}
DISCOVERY_OUTROOT=${NUMT_DISCOVERY_OUTROOT}
BESTHIT_OUTDIR=${NUMT_BESTHIT_OUTDIR}
CONFIG

    log "Starting NUMT discovery for ${SAMPLE_ID} with config ${NUMT_CONFIG}"
    (
        cd "${NUMT_DIR}"
        bash run_numt_end2end.sh --config "${NUMT_CONFIG}"
    )
fi

# Empty high-confidence NUMT BED files are valid for samples with no calls.
validate_file_exists "${NUMT_BED}"
log "NUMT validation passed for ${SAMPLE_ID}"

log "Starting round 1 for ${SAMPLE_ID}"
(
    cd "${ROUND1_LAUNCH_DIR}"
    nextflow run "${REPO_DIR}/primate_pipeline_numt_decoy_round1.nf" \
        -profile cluster \
        -resume \
        -w "${ROUND1_WORK_DIR}" \
        --sample_tsv "${SAMPLE_TSV}" \
        --outdir "${ROUND1_OUTDIR}" \
        --cram_dirs "${PRE_OUTPUT_DIR}" \
        --numt_bed_dir "${NUMT_BESTHIT_OUTDIR}"
)

ROUND1_BAM="${ROUND1_OUTDIR}/${SAMPLE_ID}/round_1/candidate_reads/${SAMPLE_ID}.with_mates.bam"
ROUND1_VCF_DIR="${ROUND1_OUTDIR}/${SAMPLE_ID}/round_1_variant_calling_decoy"
ROUND1_NUMT_FA="${ROUND1_OUTDIR}/${SAMPLE_ID}/round_1/numt_decoy_ref/${SAMPLE_ID}.original_numt.fa"
ROUND1_NUMT_VCF="${ROUND1_OUTDIR}/${SAMPLE_ID}/round_1/numt_decoy_variant_calling/${SAMPLE_ID}.numt_decoy.raw.vcf.gz"
validate_file_nonempty "${ROUND1_BAM}"
[[ -d "${ROUND1_VCF_DIR}" ]] || fail "Missing required directory: ${ROUND1_VCF_DIR}"
find "${ROUND1_VCF_DIR}" -type f -name "${SAMPLE_ID}.numt_decoy.clean.final.split.vcf" -size +0c -print -quit | grep -q . \
    || fail "Missing or empty round 1 final VCF under: ${ROUND1_VCF_DIR}"
validate_file_exists "${ROUND1_NUMT_FA}"
validate_file_nonempty "${ROUND1_NUMT_VCF}"
validate_file_nonempty "${ROUND1_NUMT_VCF}.tbi"
log "Round 1 validation passed for ${SAMPLE_ID}"

log "Starting round 2 for ${SAMPLE_ID}"
(
    cd "${ROUND2_LAUNCH_DIR}"
    nextflow run "${REPO_DIR}/primate_pipeline_round2_consensus_NUMT.nf" \
        -profile cluster \
        -resume \
        -w "${ROUND2_WORK_DIR}" \
        --sample_tsv "${SAMPLE_TSV}" \
        --outdir "${ROUND_OUTPUT_DIR}" \
        --round1_outdir "${ROUND1_OUTDIR}"
)

ROUND2_VCF_DIR="${ROUND_OUTPUT_DIR}/${SAMPLE_ID}/round_2_variant_calling_original_coords"
[[ -d "${ROUND2_VCF_DIR}" ]] || fail "Missing required directory: ${ROUND2_VCF_DIR}"
find "${ROUND2_VCF_DIR}" -type f -name "${SAMPLE_ID}.round2.original_coords.clean.final.split.vcf.gz" -size +0c -print -quit | grep -q . \
    || fail "Missing or empty round 2 final VCF under: ${ROUND2_VCF_DIR}"
log "Round 2 validation passed for ${SAMPLE_ID}"

if [[ "${CLEAN_ON_SUCCESS}" == "1" || "${CLEAN_ON_SUCCESS}" == "true" ]]; then
    log "CLEAN_ON_SUCCESS=${CLEAN_ON_SUCCESS}; deleting this sample's Nextflow work directories only"
    rm -rf "${PRE_WORK_DIR}" "${ROUND1_WORK_DIR}" "${ROUND2_WORK_DIR}"
else
    log "CLEAN_ON_SUCCESS=${CLEAN_ON_SUCCESS}; retaining work directories for resume/debugging"
fi

log "Completed per-sample streaming pipeline for ${SAMPLE_ID}"
