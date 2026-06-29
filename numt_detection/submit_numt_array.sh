#!/usr/bin/env bash
#SBATCH -J numtDisc
#SBATCH -p day
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH -o logs/numtDisc_%A_%a.out
#SBATCH -e logs/numtDisc_%A_%a.err
#SBATCH --array=1-790%50

set -euo pipefail

WORKDIR="/nfs/roberts/pi/pi_njl27/lt692/numt_discovery"
SAMPLES_TSV="/nfs/roberts/pi/pi_njl27/lt692/numt_discovery/primate_numt_test_list.txt"
RESULTS_ROOT="/nfs/roberts/pi/pi_njl27/lt692/numt_discovery/results_strict_pad"
NUCLEAR_ONLY_REF_DIR="/nfs/roberts/pi/pi_njl27/lt692/primate_ref_files/nuclear_only_refs"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

MIN_MAPQ=20
MIN_DEPTH=3
MIN_READS=5
MIN_LEN=50
MERGE_GAP=50
PAD=100

cd "${WORKDIR}"
mkdir -p logs
mkdir -p "${RESULTS_ROOT}"

module load SAMtools/1.21-GCC-13.3.0
module load BWA/0.7.18-GCCcore-13.3.0
module load miniconda/24.11.3
conda activate blast_env

echo "Job ID     : ${SLURM_JOB_ID}"
echo "Array task : ${SLURM_ARRAY_TASK_ID}"
echo "Start time : $(date)"

# 读取当前 task 对应的 sample_id（假设第1列是 sample_id）
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES_TSV}" || true)

if [[ -z "${LINE}" ]]; then
  echo "No line found for task ${SLURM_ARRAY_TASK_ID}, skip."
  exit 0
fi

SAMPLE_ID=$(printf '%s\n' "${LINE}" | cut -f1)

if [[ -z "${SAMPLE_ID}" ]]; then
  echo "Empty sample_id at task ${SLURM_ARRAY_TASK_ID}, skip."
  exit 0
fi

SAMPLE_OUTDIR="${RESULTS_ROOT}/${SAMPLE_ID}"

# 如果该样本目录下已经存在 *.numt_candidates.bed，则跳过
if compgen -G "${SAMPLE_OUTDIR}"/*.numt_candidates.bed > /dev/null; then
  echo "Skip ${SAMPLE_ID}: existing numt_candidates.bed found in ${SAMPLE_OUTDIR}"
  exit 0
fi

bash run_numt_one_from_3col.sh \
  "${SAMPLES_TSV}" \
  "${SLURM_ARRAY_TASK_ID}" \
  "${RESULTS_ROOT}" \
  "${NUCLEAR_ONLY_REF_DIR}" \
  "${THREADS}" \
  "${MIN_MAPQ}" \
  "${MIN_DEPTH}" \
  "${MIN_READS}" \
  "${MIN_LEN}" \
  "${MERGE_GAP}" \
  "${PAD}"

echo "End time   : $(date)"