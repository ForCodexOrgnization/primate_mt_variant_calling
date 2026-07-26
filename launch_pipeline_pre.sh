#!/bin/bash
#SBATCH --job-name=NF_Primate_Batch
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --requeue
#SBATCH --signal=B:USR1@300
#SBATCH --output=log_pre/nf_batch_%A_%a.log

set -euo pipefail

NEXTFLOW_PID=""
ACTIVE_BATCH_FILE=""
ACTIVE_BATCH_ID=""
REQUEUE_IN_PROGRESS=0

handle_requeue() {
    REQUEUE_IN_PROGRESS=1
    echo "INFO: Received USR1 before Slurm walltime; preserving the batch work directory."
    if [[ -n "${NEXTFLOW_PID}" ]]; then
        echo "INFO: Stopping Nextflow coordinator PID ${NEXTFLOW_PID}."
        kill -TERM "${NEXTFLOW_PID}" 2>/dev/null || true
        wait "${NEXTFLOW_PID}" 2>/dev/null || true
        NEXTFLOW_PID=""
    fi

    if [[ -n "${SLURM_JOB_ID:-}" ]] && scontrol requeue "${SLURM_JOB_ID}"; then
        echo "INFO: Requeued Slurm job ${SLURM_JOB_ID}; it will use the same batch file, work directory, and -resume."
    else
        echo "WARN: scontrol requeue was unavailable; submitting exactly one replacement job." >&2
        [[ -n "${ACTIVE_BATCH_FILE}" && -n "${ACTIVE_BATCH_ID}" ]] || {
            echo "ERROR: Cannot safely resubmit without a fixed batch file and ID." >&2
            exit 1
        }
        sbatch --export="ALL,BATCH_FILE=${ACTIVE_BATCH_FILE},BATCH_ID=${ACTIVE_BATCH_ID}" "$0"
        echo "INFO: Replacement submitted with fixed BATCH_FILE=${ACTIVE_BATCH_FILE} and BATCH_ID=${ACTIVE_BATCH_ID}."
    fi
    exit 0
}

trap handle_requeue USR1

run_nextflow() {
    set +e
    nextflow run "${SUBMIT_DIR}/preprocessing.nf" \
        -c "${NF_CONFIG_PATH}" -profile cluster -resume -w "${WORK_DIR}" \
        --sample_tsv "${ACTIVE_BATCH_FILE}" --outdir "${OUTPUT_DIR}" \
        --enable_chunked_alignment "${ENABLE_CHUNKED_ALIGNMENT:-true}" \
        --skip_existing_cram "${SKIP_EXISTING_CRAM:-true}" \
        --force_reprocess_existing_cram "${FORCE_REPROCESS_EXISTING_CRAM:-false}" \
        --samtools_bin "${SAMTOOLS_BIN:-samtools}" &
    NEXTFLOW_PID=$!
    wait "${NEXTFLOW_PID}"
    local status=$?
    NEXTFLOW_PID=""
    set -e
    return "${status}"
}

verify_batch_outputs() {
    local verify_exit=0 sample_id species ref_name cram_path crai_path marker_path
    while IFS=$'\t' read -r sample_id species ref_name _; do
        [[ -n "${sample_id}" && "${sample_id}" != "sample" && "${sample_id}" != "sample_id" ]] || continue
        cram_path="${OUTPUT_DIR}/${sample_id}/alignment/${sample_id}.cram"
        if [[ -f "${cram_path}.crai" ]]; then crai_path="${cram_path}.crai"; else crai_path="${OUTPUT_DIR}/${sample_id}/alignment/${sample_id}.crai"; fi
        marker_path="${OUTPUT_DIR}/${sample_id}/alignment/${sample_id}.cram.complete"
        if "${SUBMIT_DIR}/scripts/validate_cram.sh" --cram "${cram_path}" --crai "${crai_path}" \
             --marker "${marker_path}" --samtools "${SAMTOOLS_BIN:-samtools}" \
             --min-cram-size 1024 --min-crai-size 16; then
            :
        else
            status=$?
            echo "ERROR: CRAM validator failed for ${sample_id} (status ${status}): ${cram_path}" >&2
            verify_exit=1
        fi
    done < "${ACTIVE_BATCH_FILE}"
    return "${verify_exit}"
}

# --- 用户配置区 ---
FULL_SAMPLE_LIST="${FULL_SAMPLE_LIST:-/nfs/roberts/project/pi_njl27/lt692/primate_mt_variant_calling/human_sample.txt}"
BATCH_SIZE="${BATCH_SIZE:-5}"
CONCURRENT_BATCHES="${CONCURRENT_BATCHES:-2}"
NF_BASE_WORK_DIR="${NF_BASE_WORK_DIR:-/nfs/roberts/scratch/pi_njl27/lt692/nf_work_dir_pre}"
OUTPUT_DIR="${OUTPUT_DIR:-/nfs/roberts/scratch/pi_njl27/lt692/primate_results}"

NF_CONFIG_FILE="${NF_CONFIG_FILE:-nextflow.config}"
module load Nextflow/24.10.2
module load SAMtools/1.21-GCC-13.3.0
SAMTOOLS_BIN="${SAMTOOLS_BIN:-$(command -v samtools || true)}"
export SAMTOOLS_BIN
echo "INFO: Coordinator samtools: ${SAMTOOLS_BIN:-not found}"
if [[ -n "${SAMTOOLS_BIN}" ]]; then "${SAMTOOLS_BIN}" --version | head -n 2 || true; fi
if [[ "${SKIP_EXISTING_CRAM:-true}" == "true" && -z "${SAMTOOLS_BIN}" ]]; then
    echo "ERROR: samtools is required in the Nextflow coordinator environment to validate existing CRAM files." >&2
    exit 1
fi
# ==============================================================================

# BATCH_FILE mode: run exactly one pre-split batch directly instead of creating
# a Slurm array. This is used by launch_pipeline_all_per_batch.sh so all steps
# in one fixed batch chain use the same sample list.
if [ -n "${BATCH_FILE:-}" ]; then
    echo "--- Running in BATCH_FILE mode ---"
    [[ -s "${BATCH_FILE}" ]] || { echo "Error: BATCH_FILE does not exist or is empty: ${BATCH_FILE}" >&2; exit 1; }

    if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/preprocessing.nf" ]]; then
        SUBMIT_DIR="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
    else
        SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    fi

    RUN_DIR="${BATCH_ID:-$(basename "${BATCH_FILE}")}"
    ACTIVE_BATCH_ID="${RUN_DIR}"
    ACTIVE_BATCH_FILE="$(readlink -f "${BATCH_FILE}")"
    WORK_DIR="${NF_BASE_WORK_DIR}/${RUN_DIR}"
    mkdir -p "${WORK_DIR}"
    cd "${WORK_DIR}"

    NF_CONFIG_PATH="$(if [[ "${NF_CONFIG_FILE}" = /* ]]; then printf '%s' "${NF_CONFIG_FILE}"; else printf '%s' "${SUBMIT_DIR}/${NF_CONFIG_FILE}"; fi)"
    echo "INFO: Processing fixed batch file: ${ACTIVE_BATCH_FILE}"
    echo "INFO: Using persistent Nextflow work directory: ${WORK_DIR}"
    echo "INFO: ENABLE_CHUNKED_ALIGNMENT=${ENABLE_CHUNKED_ALIGNMENT:-true}"
    echo "INFO: Effective Nextflow config: $(if [[ "${NF_CONFIG_FILE}" = /* ]]; then printf '%s' "${NF_CONFIG_FILE}"; else printf '%s' "${SUBMIT_DIR}/${NF_CONFIG_FILE}"; fi)"

    run_nextflow || { status=$?; echo "ERROR: Nextflow failed (${status}); retaining ${WORK_DIR}." >&2; exit "${status}"; }
    verify_batch_outputs || { echo "ERROR: Output verification failed; retaining ${WORK_DIR}." >&2; exit 1; }

    echo "BATCH_FILE mode completed successfully for ${BATCH_FILE}"
    cd "${SUBMIT_DIR}"
    find "$(readlink -f "${OUTPUT_DIR}")" -type d -name "inputs" -exec rm -rf {} +
    if [[ "${DEFER_WORK_DIR_CLEANUP:-0}" == "1" ]]; then
        echo "INFO: DEFER_WORK_DIR_CLEANUP=1; retaining Nextflow work directory for parent cleanup: ${WORK_DIR}"
    else
        rm -rf "${WORK_DIR}"
    fi
    exit 0
fi

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    # --- 登录节点逻辑 (Master Mode) ---
    echo "--- Running in Master Mode on Login/Coordinator Node ---"

    # 确保 NF_BASE_WORK_DIR 存在
    mkdir -p "${NF_BASE_WORK_DIR}"

    BATCH_MANIFEST="${NF_BASE_WORK_DIR}/sample_batches.sha256"
    CURRENT_BATCH_SIGNATURE="$(sha256sum "${FULL_SAMPLE_LIST}" | awk '{print $1}'):${BATCH_SIZE}"
    if compgen -G "${NF_BASE_WORK_DIR}/sample_batch_*" >/dev/null; then
        [[ -s "${BATCH_MANIFEST}" ]] || { echo "ERROR: Existing batches lack ${BATCH_MANIFEST}; refusing an unsafe resume." >&2; exit 1; }
        [[ "$(cat "${BATCH_MANIFEST}")" == "${CURRENT_BATCH_SIGNATURE}" ]] || { echo "ERROR: Sample list or BATCH_SIZE changed; preserve the old inputs or use a new NF_BASE_WORK_DIR." >&2; exit 1; }
        echo "INFO: Reusing stable existing sample_batch_* files."
    else
        echo "Splitting ${FULL_SAMPLE_LIST} into stable batches of ${BATCH_SIZE} samples..."
        split -l "${BATCH_SIZE}" "${FULL_SAMPLE_LIST}" "${NF_BASE_WORK_DIR}/sample_batch_"
        printf '%s\n' "${CURRENT_BATCH_SIGNATURE}" > "${BATCH_MANIFEST}"
    fi

    mapfile -t BATCH_FILES < <(find "${NF_BASE_WORK_DIR}" -maxdepth 1 -type f -name 'sample_batch_*' | sort)
    NUM_BATCHES=${#BATCH_FILES[@]}
    if [ "${NUM_BATCHES}" -eq 0 ]; then
        echo "Error: No batch files were created under ${NF_BASE_WORK_DIR}."
        exit 1
    fi
    ARRAY_INDEX=$((NUM_BATCHES - 1))

    echo "Found ${NUM_BATCHES} batches. Submitting job array to run ${CONCURRENT_BATCHES} batches concurrently..."
    sbatch --array=0-${ARRAY_INDEX}%${CONCURRENT_BATCHES} "$0"

    echo "Job array submitted to Slurm. Monitor with 'squeue -u \$USER'."
    exit 0

else
    echo "================================================================="
    echo "--- Running in Worker Mode on Compute Node (Task ${SLURM_ARRAY_TASK_ID}) ---"
    echo "================================================================="

    LOG_DIR="log_pre"
    mkdir -p "$LOG_DIR"
    echo "[*] Log directory: $(pwd)/${LOG_DIR}"

    SUBMIT_DIR="$SLURM_SUBMIT_DIR"

    # 1. RUN_DIR：用于区分不同数组任务的 Nextflow run
    RUN_DIR="batch_${SLURM_ARRAY_TASK_ID}"

    # 2. WORK_DIR：每个数组任务自己的 Nextflow work 目录（固定且唯一）
    WORK_DIR="${NF_BASE_WORK_DIR}/${RUN_DIR}"

    # 确保 Nextflow work 目录存在
    mkdir -p "${WORK_DIR}"
    cd "${WORK_DIR}"

    TASK_ID_PLUS_ONE=$((SLURM_ARRAY_TASK_ID + 1))

    # 从 NF_BASE_WORK_DIR 下的 sample_batch_* 中，按排序顺序取第 N 个
    mapfile -t BATCH_FILES < <(find "${NF_BASE_WORK_DIR}" -maxdepth 1 -type f -name 'sample_batch_*' | sort)
    BATCH_FILE="${BATCH_FILES[${SLURM_ARRAY_TASK_ID}]:-}"

    if [ -z "${BATCH_FILE}" ]; then
        echo "Error: Could not find batch file for task ID ${SLURM_ARRAY_TASK_ID} in ${NF_BASE_WORK_DIR}."
        exit 1
    fi
    ACTIVE_BATCH_FILE="$(readlink -f "${BATCH_FILE}")"
    ACTIVE_BATCH_ID="${RUN_DIR}"
    NF_CONFIG_PATH="$(if [[ "${NF_CONFIG_FILE}" = /* ]]; then printf '%s' "${NF_CONFIG_FILE}"; else printf '%s' "${SUBMIT_DIR}/${NF_CONFIG_FILE}"; fi)"

    echo "INFO: This task will process batch file: ${BATCH_FILE}"
    echo "INFO: Using persistent Nextflow work directory: ${WORK_DIR}"
    echo "INFO: ENABLE_CHUNKED_ALIGNMENT=${ENABLE_CHUNKED_ALIGNMENT:-true}"
    echo "INFO: Effective Nextflow config: $(if [[ "${NF_CONFIG_FILE}" = /* ]]; then printf '%s' "${NF_CONFIG_FILE}"; else printf '%s' "${SUBMIT_DIR}/${NF_CONFIG_FILE}"; fi)"

    if run_nextflow; then NF_EXIT=0; else NF_EXIT=$?; fi

    if [ ${NF_EXIT} -eq 0 ]; then
        echo "Verifying CRAM/CRAI outputs for batch ${SLURM_ARRAY_TASK_ID} before cleanup..."
        if ! verify_batch_outputs; then
            echo "Batch ${SLURM_ARRAY_TASK_ID} failed output verification; retaining work directory for debugging: ${WORK_DIR}" >&2
            exit 1
        fi

        echo "Batch ${SLURM_ARRAY_TASK_ID} completed successfully."
        # 成功：清理 Nextflow work 目录以释放空间
        echo "Cleaning up Nextflow work directory: ${WORK_DIR}"
        cd "${SUBMIT_DIR}"
        
        # 自定义清理逻辑：删除输出目录里的 inputs 目录
        find "$(readlink -f "${OUTPUT_DIR}")" -type d -name "inputs" -exec rm -rf {} +
        rm -rf "${WORK_DIR}"
        #rm -f "${BATCH_FILE}"

    else
        echo "Batch ${SLURM_ARRAY_TASK_ID} failed (exit code ${NF_EXIT})."
        # 失败：保留 work 目录，以便下次 -resume / debug
        echo "Retaining work directory for next run/debugging: ${WORK_DIR}"
    fi

    echo "--- Finished Job Array Task ${SLURM_ARRAY_TASK_ID} ---"
fi
