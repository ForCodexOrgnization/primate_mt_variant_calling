#!/usr/bin/env bash
#SBATCH --job-name=watch_pre_done
#SBATCH --output=log_streaming/watch_pre_done_%j.log
#SBATCH --error=log_streaming/watch_pre_done_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -p ycga

set -euo pipefail

module load SAMtools/1.21-GCC-12.2.0 || true

PRE_OUTPUT_DIR="/home/lt692/scratch_pi_njl27/lt692/primate_results_test" \
ROUND_OUTPUT_DIR="/home/lt692/scratch_pi_njl27/lt692/primate_results_test" \
MASTER_SAMPLE_LIST="/home/lt692/project_pi_njl27/lt692/primate_mt_variant_calling/human_sample.txt" \
REPO_DIR="/home/lt692/project_pi_njl27/lt692/primate_mt_variant_calling" \
NF_BASE_WORK_DIR="/home/lt692/ycga_work/nf_work_dir_streaming_per_sample" \
CLEAN_ON_SUCCESS=1 \
MAX_CONCURRENT=1 \
WATCH_INTERVAL_SECONDS=300 \
NF_CONFIG_FILE="${NF_CONFIG_FILE:-nextflow_mcc.config}" \
bash /home/lt692/project_pi_njl27/lt692/primate_mt_variant_calling/watch_completed_pre_and_launch_rounds.sh
