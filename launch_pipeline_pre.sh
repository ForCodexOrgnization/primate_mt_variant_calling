#!/bin/bash
#SBATCH --job-name=NF_Primate_Batch
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=log_pre/nf_batch_%A_%a.log

set -e

# --- 用户配置区 ---
FULL_SAMPLE_LIST="${FULL_SAMPLE_LIST:-/nfs/roberts/project/pi_njl27/lt692/primate_mt_variant_calling/human_sample.txt}"
BATCH_SIZE="${BATCH_SIZE:-5}"
CONCURRENT_BATCHES="${CONCURRENT_BATCHES:-2}"
NF_BASE_WORK_DIR="${NF_BASE_WORK_DIR:-/nfs/roberts/scratch/pi_njl27/lt692/nf_work_dir_pre}"
OUTPUT_DIR="${OUTPUT_DIR:-/nfs/roberts/scratch/pi_njl27/lt692/primate_results}"

module load Nextflow/24.10.2
# ==============================================================================

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    # --- 登录节点逻辑 (Master Mode) ---
    echo "--- Running in Master Mode on Login/Coordinator Node ---"

    # 确保 NF_BASE_WORK_DIR 存在
    mkdir -p "${NF_BASE_WORK_DIR}"

    echo "Cleaning up old batch files in ${NF_BASE_WORK_DIR}."
    # 只清 sample_batch_*，不动 batch_* 目录，方便 Nextflow -resume
    rm -f "${NF_BASE_WORK_DIR}/sample_batch_"*

    echo "Splitting ${FULL_SAMPLE_LIST} into batches of ${BATCH_SIZE} samples..."
    # 直接把 batch 文件写到 NF_BASE_WORK_DIR 下
    split -l "${BATCH_SIZE}" "${FULL_SAMPLE_LIST}" "${NF_BASE_WORK_DIR}/sample_batch_"

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

    echo "INFO: This task will process batch file: ${BATCH_FILE}"
    echo "INFO: Using persistent Nextflow work directory: ${WORK_DIR}"

    set +e
    nextflow run "${SUBMIT_DIR}/preprocessing.nf" \
        -profile cluster \
        -resume \
        -w "${WORK_DIR}" \
        --sample_tsv "${BATCH_FILE}" \
        --outdir "${OUTPUT_DIR}"

    NF_EXIT=$?
    set -e

    if [ ${NF_EXIT} -eq 0 ]; then
        echo "Verifying CRAM/CRAI outputs for batch ${SLURM_ARRAY_TASK_ID} before cleanup..."
        VERIFY_EXIT=0
        while IFS=$'\t' read -r SAMPLE_ID _; do
            # Skip blank lines and an optional header.
            if [ -z "${SAMPLE_ID}" ] || [ "${SAMPLE_ID}" = "sample" ] || [ "${SAMPLE_ID}" = "sample_id" ]; then
                continue
            fi
            CRAM_PATH="${OUTPUT_DIR}/${SAMPLE_ID}/alignment/${SAMPLE_ID}.cram"
            CRAI_PATH="${OUTPUT_DIR}/${SAMPLE_ID}/alignment/${SAMPLE_ID}.cram.crai"
            if [ ! -s "${CRAM_PATH}" ] || [ ! -s "${CRAI_PATH}" ]; then
                echo "ERROR: Nextflow exited 0 but expected CRAM/CRAI is missing or empty for sample '${SAMPLE_ID}'." >&2
                echo "       CRAM: ${CRAM_PATH}" >&2
                echo "       CRAI: ${CRAI_PATH}" >&2
                VERIFY_EXIT=1
            fi
        done < "${BATCH_FILE}"

        if [ ${VERIFY_EXIT} -ne 0 ]; then
            echo "Batch ${SLURM_ARRAY_TASK_ID} failed output verification; retaining work directory for debugging: ${WORK_DIR}" >&2
            exit ${VERIFY_EXIT}
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
