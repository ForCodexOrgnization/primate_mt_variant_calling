#!/bin/bash
#SBATCH --job-name=NF_Primate_Batch
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=log_round2/nf_batch_%A_%a.log

set -euo pipefail

# --- 用户配置区 ---
FULL_SAMPLE_LIST="${FULL_SAMPLE_LIST:-/home/lt692/project_pi_njl27/lt692/primate_mt_variant_calling/human_sample.txt}"
BATCH_SIZE="${BATCH_SIZE:-5}"
CONCURRENT_BATCHES="${CONCURRENT_BATCHES:-2}"
NF_BASE_WORK_DIR="${NF_BASE_WORK_DIR:-/nfs/roberts/scratch/pi_njl27/lt692/nf_work_dir_round2}"
if [[ -n "${BATCH_FILE:-}" ]]; then
    OUTPUT_DIR="${OUTPUT_DIR:?ERROR: OUTPUT_DIR must be explicitly set in BATCH_FILE mode}"
    ROUND1_OUTDIR="${ROUND1_OUTDIR:?ERROR: ROUND1_OUTDIR must be explicitly set in BATCH_FILE mode}"
else
    OUTPUT_DIR="${OUTPUT_DIR:-/nfs/roberts/scratch/pi_njl27/lt692/primate_results}"
    ROUND1_OUTDIR="${ROUND1_OUTDIR:-${OUTPUT_DIR}}"
fi

NF_CONFIG_FILE="${NF_CONFIG_FILE:-nextflow.config}"
module load Nextflow/24.10.2
# ==============================================================================

# BATCH_FILE mode: run exactly one pre-split batch directly instead of creating
# a Slurm array. This is used by launch_pipeline_all_per_batch.sh so all steps
# in one fixed batch chain use the same sample list.
if [ -n "${BATCH_FILE:-}" ]; then
    echo "--- Running in BATCH_FILE mode ---"
    [[ -s "${BATCH_FILE}" ]] || { echo "Error: BATCH_FILE does not exist or is empty: ${BATCH_FILE}" >&2; exit 1; }

    if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/primate_pipeline_round2_consensus_NUMT.nf" ]]; then
        SUBMIT_DIR="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
    else
        SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    fi

    RUN_DIR="${BATCH_ID:-$(basename "${BATCH_FILE}")}"
    WORK_DIR="${NF_BASE_WORK_DIR}/${RUN_DIR}"
    mkdir -p "${WORK_DIR}"
    cd "${WORK_DIR}"

    echo "INFO: Processing fixed batch file: ${BATCH_FILE}"
    echo "INFO: Using persistent Nextflow work directory: ${WORK_DIR}"
    echo "INFO: round2 BATCH_FILE mode parameters:"
    echo "  BATCH_FILE=${BATCH_FILE}"
    echo "  OUTPUT_DIR=${OUTPUT_DIR}"
    echo "  ROUND1_OUTDIR=${ROUND1_OUTDIR}"
    echo "  NF_BASE_WORK_DIR=${NF_BASE_WORK_DIR}"
    echo "  WORK_DIR=${WORK_DIR}"

    nextflow run "${SUBMIT_DIR}/primate_pipeline_round2_consensus_NUMT.nf" \
        -c "$(if [[ "${NF_CONFIG_FILE}" = /* ]]; then printf '%s' "${NF_CONFIG_FILE}"; else printf '%s' "${SUBMIT_DIR}/${NF_CONFIG_FILE}"; fi)" \
        -profile cluster \
        -resume \
        -w "${WORK_DIR}" \
        --sample_tsv "${BATCH_FILE}" \
        --outdir "${OUTPUT_DIR}" \
        --round1_outdir "${ROUND1_OUTDIR}"

    echo "BATCH_FILE mode completed successfully for ${BATCH_FILE}"
    cd "${SUBMIT_DIR}"
    find "$(readlink -f "${OUTPUT_DIR}")" -type d -name "inputs" -exec rm -rf {} +
    rm -rf "${WORK_DIR}"
    exit 0
fi

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "--- Running in Master Mode on Login/Coordinator Node ---"

    mkdir -p "${NF_BASE_WORK_DIR}"
    mkdir -p log_round2

    echo "Cleaning up old batch files in ${NF_BASE_WORK_DIR}"
    rm -f "${NF_BASE_WORK_DIR}/sample_batch_"*

    echo "Splitting ${FULL_SAMPLE_LIST} into batches of ${BATCH_SIZE} samples..."
    split -l "${BATCH_SIZE}" "${FULL_SAMPLE_LIST}" "${NF_BASE_WORK_DIR}/sample_batch_"

    mapfile -t BATCH_FILES < <(find "${NF_BASE_WORK_DIR}" -maxdepth 1 -type f -name "sample_batch_*" | sort)
    NUM_BATCHES=${#BATCH_FILES[@]}
    if [ "${NUM_BATCHES}" -eq 0 ]; then
        echo "Error: No batch files were created under ${NF_BASE_WORK_DIR}"
        exit 1
    fi

    ARRAY_INDEX=$((NUM_BATCHES - 1))

    echo "Found ${NUM_BATCHES} batches. Submitting job array with concurrency ${CONCURRENT_BATCHES}..."
    sbatch --array=0-"${ARRAY_INDEX}"%"${CONCURRENT_BATCHES}" "$0"

    echo "Job array submitted. Monitor with: squeue -u \$USER"
    exit 0

else
    echo "================================================================="
    echo "--- Running in Worker Mode on Compute Node (Task ${SLURM_ARRAY_TASK_ID}) ---"
    echo "================================================================="

    LOG_DIR="${SLURM_SUBMIT_DIR}/log_round2"
    mkdir -p "${LOG_DIR}"

    SUBMIT_DIR="${SLURM_SUBMIT_DIR}"
    RUN_DIR="batch_${SLURM_ARRAY_TASK_ID}"
    WORK_DIR="${NF_BASE_WORK_DIR}/${RUN_DIR}"

    mkdir -p "${WORK_DIR}"
    cd "${WORK_DIR}"

    mapfile -t BATCH_FILES < <(find "${NF_BASE_WORK_DIR}" -maxdepth 1 -type f -name "sample_batch_*" | sort)
    BATCH_FILE="${BATCH_FILES[${SLURM_ARRAY_TASK_ID}]:-}"

    if [ -z "${BATCH_FILE}" ]; then
        echo "Error: Could not find batch file for task ID ${SLURM_ARRAY_TASK_ID} in ${NF_BASE_WORK_DIR}"
        exit 1
    fi

    echo "INFO: This task will process batch file: ${BATCH_FILE}"
    echo "INFO: Using persistent Nextflow work directory: ${WORK_DIR}"

    set +e
    nextflow run "${SUBMIT_DIR}/primate_pipeline_round2_consensus_NUMT.nf" \
        -c "$(if [[ "${NF_CONFIG_FILE}" = /* ]]; then printf '%s' "${NF_CONFIG_FILE}"; else printf '%s' "${SUBMIT_DIR}/${NF_CONFIG_FILE}"; fi)" \
        -profile cluster \
        -resume \
        -w "${WORK_DIR}" \
        --sample_tsv "${BATCH_FILE}" \
        --outdir "${OUTPUT_DIR}" \
        --round1_outdir "${ROUND1_OUTDIR}"

    NF_EXIT=$?
    set -e

    if [ "${NF_EXIT}" -eq 0 ]; then
        echo "Batch ${SLURM_ARRAY_TASK_ID} completed successfully."
        echo "Cleaning up Nextflow work directory: ${WORK_DIR}"

        cd "${SUBMIT_DIR}"

        # 谨慎保留；确认你确实想删输出目录下所有 inputs 目录
        find "$(readlink -f "${OUTPUT_DIR}")" -type d -name "inputs" -exec rm -rf {} +

        rm -rf "${WORK_DIR}"
        # rm -f "${BATCH_FILE}"
    else
        echo "Batch ${SLURM_ARRAY_TASK_ID} failed (exit code ${NF_EXIT})."
        echo "Retaining work directory for debugging/resume: ${WORK_DIR}"
        exit "${NF_EXIT}"
    fi

    echo "--- Finished Job Array Task ${SLURM_ARRAY_TASK_ID} ---"
fi
