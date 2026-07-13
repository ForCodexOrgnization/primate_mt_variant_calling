#!/usr/bin/env bash
#SBATCH --job-name=NF_Primate_Chain
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=log_all_per_batch/nf_chain_%A_%a.log

set -euo pipefail

# ==============================================================================
# Per-batch chained launcher: each fixed batch runs pre -> NUMT -> round1 -> round2
# without waiting for other batches to finish the same step.
#
# Example:
#   FULL_SAMPLE_LIST=/path/samples.tsv BATCH_SIZE=5 CHAIN_CONCURRENT_BATCHES=2 \
#     bash launch_pipeline_all_per_batch.sh
# ==============================================================================

FULL_SAMPLE_LIST="${FULL_SAMPLE_LIST:-/nfs/roberts/project/pi_njl27/lt692/primate_mt_variant_calling/bouchet_sample_list.txt}"
PRE_OUTPUT_DIR="${PRE_OUTPUT_DIR:-/nfs/roberts/pi/pi_njl27/lt692/primate_results}"
ROUND_OUTPUT_DIR="${ROUND_OUTPUT_DIR:-/nfs/roberts/pi/pi_njl27/lt692/primate_results}"
NF_BASE_WORK_DIR="${NF_BASE_WORK_DIR:-/home/lt692/scratch_pi_njl27/lt692/nf_work_dir_all_per_batch}"
PRE_NF_BASE_WORK_DIR="${PRE_NF_BASE_WORK_DIR:-${NF_BASE_WORK_DIR}/pre}"
ROUND1_NF_BASE_WORK_DIR="${ROUND1_NF_BASE_WORK_DIR:-${NF_BASE_WORK_DIR}/round1}"
ROUND2_NF_BASE_WORK_DIR="${ROUND2_NF_BASE_WORK_DIR:-${NF_BASE_WORK_DIR}/round2}"
NUMT_NF_BASE_WORK_DIR="${NUMT_NF_BASE_WORK_DIR:-${NF_BASE_WORK_DIR}/numt}"
ROUND1_OUTDIR="${ROUND1_OUTDIR:-${ROUND_OUTPUT_DIR}}"

PRE_LAUNCH_SCRIPT="${PRE_LAUNCH_SCRIPT:-launch_pipeline_pre.sh}"
NUMT_LAUNCH_SCRIPT="${NUMT_LAUNCH_SCRIPT:-numt_detection/numt_end2end.nf}"
ROUND1_LAUNCH_SCRIPT="${ROUND1_LAUNCH_SCRIPT:-launch_pipeline_round1.sh}"
ROUND2_LAUNCH_SCRIPT="${ROUND2_LAUNCH_SCRIPT:-launch_pipeline_round2.sh}"

NUMT_DISCOVERY_OUTROOT="${NUMT_DISCOVERY_OUTROOT:-${ROUND_OUTPUT_DIR}/numt_discovery}"
NUMT_BESTHIT_OUTDIR="${NUMT_BESTHIT_OUTDIR:-${ROUND_OUTPUT_DIR}/numt_besthit}"
REF_DIR="${REF_DIR:-/home/lt692/scratch_pi_njl27/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_chrM}"
GLOBAL_REF_DIR="${GLOBAL_REF_DIR:-/home/lt692/scratch_pi_njl27/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_whole}"
NUCLEAR_ONLY_REF_DIR="${NUCLEAR_ONLY_REF_DIR:-/home/lt692/scratch_pi_njl27/lt692/primate_mtDNA_analysis/references/variant_calling/nuclear_only_refs}"
BATCH_SIZE="${BATCH_SIZE:-10}"
CHAIN_CONCURRENT_BATCHES="${CHAIN_CONCURRENT_BATCHES:-2}"
NUMT_CONCURRENT="${NUMT_CONCURRENT:-${CONCURRENT:-${CHAIN_CONCURRENT_BATCHES}}}"
POLL_SECONDS="${POLL_SECONDS:-120}"
LOG_DIR="${LOG_DIR:-log_all_per_batch}"
BATCH_LIST_DIR="${BATCH_LIST_DIR:-${NF_BASE_WORK_DIR}/batch_lists}"
NF_CONFIG_FILE="${NF_CONFIG_FILE:-nextflow.config}"
CLEAN_ON_SUCCESS="${CLEAN_ON_SUCCESS:-1}"

log() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" >&2; }
fail() { echo "ERROR: $*" >&2; exit 1; }

if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/preprocessing.nf" ]]; then
    repo_dir="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
else
    repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
cd "${repo_dir}"
NF_CONFIG_PATH="$(if [[ "${NF_CONFIG_FILE}" = /* ]]; then printf '%s' "${NF_CONFIG_FILE}"; else printf '%s' "${repo_dir}/${NF_CONFIG_FILE}"; fi)"
[[ -s "${NF_CONFIG_PATH}" ]] || fail "NF_CONFIG_FILE does not exist or is empty: ${NF_CONFIG_PATH}"
mkdir -p "${LOG_DIR}"

echo "INFO: FULL_SAMPLE_LIST=${FULL_SAMPLE_LIST}"
echo "INFO: PRE_OUTPUT_DIR=${PRE_OUTPUT_DIR}"
echo "INFO: ROUND_OUTPUT_DIR=${ROUND_OUTPUT_DIR}"
echo "INFO: ROUND1_OUTDIR=${ROUND1_OUTDIR}"
echo "INFO: NUMT_DISCOVERY_OUTROOT=${NUMT_DISCOVERY_OUTROOT}"
echo "INFO: NUMT_BESTHIT_OUTDIR=${NUMT_BESTHIT_OUTDIR}"
echo "INFO: NUMT_NF_BASE_WORK_DIR=${NUMT_NF_BASE_WORK_DIR}"
echo "INFO: NF_BASE_WORK_DIR=${NF_BASE_WORK_DIR}"
echo "INFO: NF_CONFIG_FILE=${NF_CONFIG_FILE}"
echo "INFO: NF_CONFIG_PATH=${NF_CONFIG_PATH}"
echo "INFO: CLEAN_ON_SUCCESS=${CLEAN_ON_SUCCESS}"

cleanup_work_dir_if_requested() {
    local stage_name="$1"
    local work_dir="$2"

    if [[ "${CLEAN_ON_SUCCESS}" == "1" || "${CLEAN_ON_SUCCESS}" == "true" ]]; then
        if [[ -n "${work_dir}" && -d "${work_dir}" ]]; then
            log "CLEAN_ON_SUCCESS=${CLEAN_ON_SUCCESS}; deleting ${stage_name} work directory: ${work_dir}"
            rm -rf "${work_dir}"
        else
            log "CLEAN_ON_SUCCESS=${CLEAN_ON_SUCCESS}; ${stage_name} work directory does not exist or was already removed: ${work_dir}"
        fi
    else
        log "CLEAN_ON_SUCCESS=${CLEAN_ON_SUCCESS}; retaining ${stage_name} work directory for resume/debugging: ${work_dir}"
    fi
}

read_sample_ids() {
    local sample_file="$1"
    awk -F '\t' 'NF >= 2 && $1 !~ /^#/ && $1 != "" && tolower($1) != "sample" && tolower($1) != "sample_id" {print $1}' "${sample_file}"
}

job_state_sacct() {
    local job_id="$1"
    sacct -j "${job_id}" --format=State --noheader 2>/dev/null | awk 'NF {print $1}' | sed 's/+$//' | sort -u | paste -sd, -
}

wait_for_job() {
    local job_id="$1" step_name="$2"
    log "Waiting for ${step_name} job ${job_id}"
    while true; do
        if command -v sacct >/dev/null 2>&1; then
            local states
            states="$(job_state_sacct "${job_id}")"
            if [[ -n "${states}" ]]; then
                if [[ "${states}" == *FAILED* || "${states}" == *CANCELLED* || "${states}" == *TIMEOUT* || "${states}" == *OUT_OF_MEMORY* || "${states}" == *NODE_FAIL* || "${states}" == *PREEMPTED* ]]; then
                    fail "${step_name} job ${job_id} finished unsuccessfully with states: ${states}"
                fi
                if [[ "${states}" == *COMPLETED* && "${states}" != *RUNNING* && "${states}" != *PENDING* && "${states}" != *CONFIGURING* && "${states}" != *COMPLETING* ]]; then
                    log "${step_name} job ${job_id} completed with states: ${states}"
                    return 0
                fi
                log "${step_name} job ${job_id} states: ${states}; sleeping ${POLL_SECONDS}s"
            fi
        elif ! squeue -h -j "${job_id}" >/dev/null 2>&1 || [[ -z "$(squeue -h -j "${job_id}" 2>/dev/null)" ]]; then
            log "${step_name} job ${job_id} is no longer in squeue; sacct unavailable, assuming completion"
            return 0
        else
            log "${step_name} job ${job_id} still present in squeue; sleeping ${POLL_SECONDS}s"
        fi
        sleep "${POLL_SECONDS}"
    done
}

write_numt_config() {
    local batch_file="$1" batch_id="$2" clean_samples="$3" auto_config="$4"

    awk -F '	' '
BEGIN { OFS="	" }
{
  gsub(/\r/, "", $0)
}
$0 ~ /^[[:space:]]*$/ { next }
$1 ~ /^#/ { next }
tolower($1) == "sample" { next }
tolower($1) == "sample_id" { next }

NF >= 3 {
  print $1, $2, $3
  next
}

NF == 2 {
  print $1, $2, $2
  next
}

{
  print "ERROR: malformed sample row: " $0 > "/dev/stderr"
  bad=1
}

END {
  if (bad) exit 1
}
' "${batch_file}" > "${clean_samples}"

    [[ -s "${clean_samples}" ]] || fail "No valid sample rows after cleaning ${batch_file}"

    cat > "${auto_config}" <<CONFIG
# Auto-generated by primate_mt_variant_calling/launch_pipeline_all_per_batch.sh for batch_${batch_id}
SAMPLES_TSV=${clean_samples}
CRAM_ROOT_1=${PRE_OUTPUT_DIR}
CRAM_ROOT_2=${PRE_OUTPUT_DIR}
WHOLE_REF_DIR=${GLOBAL_REF_DIR}
NUCLEAR_ONLY_REF_DIR=${NUCLEAR_ONLY_REF_DIR}
CHRM_REF_DIR=${REF_DIR}
DISCOVERY_OUTROOT=${NUMT_DISCOVERY_OUTROOT}
BESTHIT_OUTDIR=${NUMT_BESTHIT_OUTDIR}
CONFIG
}

run_numt_nextflow() {
    local batch_file="$1" batch_id="$2" batch_name="batch_${batch_id}"
    local numt_batch_work_dir="${NUMT_NF_BASE_WORK_DIR}/${batch_name}"
    local numt_batch_dir="${NUMT_DISCOVERY_OUTROOT}/nextflow_${batch_name}"
    local clean_samples="${numt_batch_dir}/numt_samples.${batch_name}.tsv"
    local auto_config="${numt_batch_dir}/numt_pipeline.${batch_name}.config"

    mkdir -p "${numt_batch_dir}" "${numt_batch_work_dir}" "${NUMT_DISCOVERY_OUTROOT}" "${NUMT_BESTHIT_OUTDIR}"
    write_numt_config "${batch_file}" "${batch_id}" "${clean_samples}" "${auto_config}"

    log "Launching NUMT Nextflow for ${batch_id}"
    log "  NUMT workflow=${NUMT_LAUNCH_SCRIPT}"
    log "  NUMT config=${auto_config}"
    log "  NUMT work dir=${numt_batch_work_dir}"
    log "  NUMT concurrency/queue size=${NUMT_CONCURRENT}"

    module load Nextflow/24.10.2 || module load Nextflow/24.04.4
    command -v nextflow >/dev/null 2>&1 || fail "nextflow not found after loading Nextflow module"
    log "  nextflow=$(command -v nextflow)"
    nextflow -version

    NXF_OPTS="${NXF_OPTS:-}" nextflow -C "${NF_CONFIG_PATH}" run "${NUMT_LAUNCH_SCRIPT}" \
        -profile cluster \
        -work-dir "${numt_batch_work_dir}" \
        -with-report "${numt_batch_dir}/nextflow.report.html" \
        -with-timeline "${numt_batch_dir}/nextflow.timeline.html" \
        -with-trace "${numt_batch_dir}/nextflow.trace.tsv" \
        -queue-size "${NUMT_CONCURRENT}" \
        --numt_config "${auto_config}"
}

validate_pre_to_round1() {
    local batch_file="$1" missing=0 sample cram crai
    log "Validating preprocessing outputs for ${batch_file}"
    while IFS= read -r sample; do
        cram="${PRE_OUTPUT_DIR}/${sample}/alignment/${sample}.cram"; crai="${cram}.crai"
        [[ -s "${cram}" ]] || { echo "MISSING/EMPTY: ${cram}" >&2; missing=1; }
        [[ -s "${crai}" ]] || { echo "MISSING/EMPTY: ${crai}" >&2; missing=1; }
    done < <(read_sample_ids "${batch_file}")
    [[ "${missing}" -eq 0 ]] || fail "Preprocessing outputs are incomplete for ${batch_file}"
}

validate_numt_to_round1() {
    local batch_file="$1" missing=0 sample bed
    log "Validating NUMT outputs for ${batch_file}"
    while IFS= read -r sample; do
        bed="${NUMT_BESTHIT_OUTDIR}/${sample}.highconf_numt.bed"
        [[ -e "${bed}" ]] || { echo "MISSING: ${bed}" >&2; missing=1; }
    done < <(read_sample_ids "${batch_file}")
    [[ "${missing}" -eq 0 ]] || fail "NUMT outputs are incomplete for ${batch_file}"
}

validate_round1_to_round2() {
    local batch_file="$1" missing=0 sample
    log "Validating round 1 outputs for ${batch_file}"
    while IFS= read -r sample; do
        local bam="${ROUND1_OUTDIR}/${sample}/round_1/candidate_reads/${sample}.with_mates.bam"
        local vcf_dir="${ROUND1_OUTDIR}/${sample}/round_1_variant_calling_decoy"
        local numt_fa="${ROUND1_OUTDIR}/${sample}/round_1/numt_decoy_ref/${sample}.original_numt.fa"
        local numt_vcf="${ROUND1_OUTDIR}/${sample}/round_1/numt_decoy_variant_calling/${sample}.numt_decoy.raw.vcf.gz"
        [[ -s "${bam}" ]] || { echo "MISSING/EMPTY: ${bam}" >&2; missing=1; }
        [[ -d "${vcf_dir}" ]] || { echo "MISSING DIR: ${vcf_dir}" >&2; missing=1; }
        find "${vcf_dir}" -type f -name "${sample}.numt_decoy.clean.final.split.vcf" -size +0c -print -quit 2>/dev/null | grep -q . || { echo "MISSING/EMPTY round1 final VCF under: ${vcf_dir}" >&2; missing=1; }
        [[ -e "${numt_fa}" ]] || { echo "MISSING: ${numt_fa}" >&2; missing=1; }
        [[ -s "${numt_vcf}" ]] || { echo "MISSING/EMPTY: ${numt_vcf}" >&2; missing=1; }
        [[ -s "${numt_vcf}.tbi" ]] || { echo "MISSING/EMPTY: ${numt_vcf}.tbi" >&2; missing=1; }
    done < <(read_sample_ids "${batch_file}")
    [[ "${missing}" -eq 0 ]] || fail "Round 1 outputs are incomplete for ${batch_file}"
}

validate_round2_final() {
    local batch_file="$1" missing=0 sample out_dir
    log "Validating round 2 outputs for ${batch_file}"
    while IFS= read -r sample; do
        out_dir="${ROUND_OUTPUT_DIR}/${sample}/round_2_variant_calling_original_coords"
        [[ -d "${out_dir}" ]] || { echo "MISSING DIR: ${out_dir}" >&2; missing=1; continue; }
        find "${out_dir}" -type f -name "${sample}.round2.original_coords.clean.final.split.vcf.gz" -size +0c -print -quit 2>/dev/null | grep -q . || { echo "MISSING/EMPTY round2 final VCF under: ${out_dir}" >&2; missing=1; }
    done < <(read_sample_ids "${batch_file}")
    [[ "${missing}" -eq 0 ]] || fail "Round 2 outputs are incomplete for ${batch_file}"
}

run_chain() {
    local batch_file="$1" batch_id="$2"
    local batch_name="batch_${batch_id}"
    local pre_batch_work_dir="${PRE_NF_BASE_WORK_DIR}/${batch_name}"
    local round1_batch_work_dir="${ROUND1_NF_BASE_WORK_DIR}/${batch_name}"
    local round2_batch_work_dir="${ROUND2_NF_BASE_WORK_DIR}/${batch_name}"

    log "Starting per-batch chain ${batch_id}: ${batch_file}"
    env BATCH_FILE="${batch_file}" BATCH_ID="${batch_name}" FULL_SAMPLE_LIST="${batch_file}" OUTPUT_DIR="${PRE_OUTPUT_DIR}" NF_BASE_WORK_DIR="${PRE_NF_BASE_WORK_DIR}" NF_CONFIG_FILE="${NF_CONFIG_FILE}" DEFER_WORK_DIR_CLEANUP=1 bash "${PRE_LAUNCH_SCRIPT}"
    validate_pre_to_round1 "${batch_file}"
    cleanup_work_dir_if_requested "pre_${batch_name}" "${pre_batch_work_dir}"
    run_numt_nextflow "${batch_file}" "${batch_id}"
    cleanup_work_dir_if_requested "numt_${batch_name}" "${NUMT_NF_BASE_WORK_DIR}/${batch_name}"
    validate_numt_to_round1 "${batch_file}"
    log "Launching round1 for ${batch_id}"
    log "  BATCH_FILE=${batch_file}"
    log "  OUTPUT_DIR=${ROUND1_OUTDIR}"
    log "  CRAM_DIRS=${PRE_OUTPUT_DIR}"
    log "  NUMT_BED_DIR=${NUMT_BESTHIT_OUTDIR}"
    log "  NF_BASE_WORK_DIR=${ROUND1_NF_BASE_WORK_DIR}"
    env BATCH_FILE="${batch_file}" BATCH_ID="${batch_name}" FULL_SAMPLE_LIST="${batch_file}" OUTPUT_DIR="${ROUND1_OUTDIR}" CRAM_DIRS="${PRE_OUTPUT_DIR}" NUMT_BED_DIR="${NUMT_BESTHIT_OUTDIR}" NF_BASE_WORK_DIR="${ROUND1_NF_BASE_WORK_DIR}" NF_CONFIG_FILE="${NF_CONFIG_FILE}" DEFER_WORK_DIR_CLEANUP=1 bash "${ROUND1_LAUNCH_SCRIPT}"
    validate_round1_to_round2 "${batch_file}"
    cleanup_work_dir_if_requested "round1_${batch_name}" "${round1_batch_work_dir}"
    log "Launching round2 for ${batch_id}"
    log "  BATCH_FILE=${batch_file}"
    log "  OUTPUT_DIR=${ROUND_OUTPUT_DIR}"
    log "  ROUND1_OUTDIR=${ROUND1_OUTDIR}"
    log "  NF_BASE_WORK_DIR=${ROUND2_NF_BASE_WORK_DIR}"
    env BATCH_FILE="${batch_file}" BATCH_ID="${batch_name}" FULL_SAMPLE_LIST="${batch_file}" OUTPUT_DIR="${ROUND_OUTPUT_DIR}" ROUND1_OUTDIR="${ROUND1_OUTDIR}" NF_BASE_WORK_DIR="${ROUND2_NF_BASE_WORK_DIR}" NF_CONFIG_FILE="${NF_CONFIG_FILE}" DEFER_WORK_DIR_CLEANUP=1 bash "${ROUND2_LAUNCH_SCRIPT}"
    validate_round2_final "${batch_file}"
    cleanup_work_dir_if_requested "round2_${batch_name}" "${round2_batch_work_dir}"
    log "Completed per-batch chain ${batch_id}"
}

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    [[ -s "${FULL_SAMPLE_LIST}" ]] || fail "FULL_SAMPLE_LIST does not exist or is empty: ${FULL_SAMPLE_LIST}"
    mkdir -p "${BATCH_LIST_DIR}" "${PRE_NF_BASE_WORK_DIR}" "${NUMT_NF_BASE_WORK_DIR}" "${ROUND1_NF_BASE_WORK_DIR}" "${ROUND2_NF_BASE_WORK_DIR}"
    rm -f "${BATCH_LIST_DIR}/sample_batch_"*
    split -l "${BATCH_SIZE}" "${FULL_SAMPLE_LIST}" "${BATCH_LIST_DIR}/sample_batch_"
    mapfile -t BATCH_FILES < <(find "${BATCH_LIST_DIR}" -maxdepth 1 -type f -name 'sample_batch_*' | sort)
    NUM_BATCHES=${#BATCH_FILES[@]}
    [[ "${NUM_BATCHES}" -gt 0 ]] || fail "No batch files were created under ${BATCH_LIST_DIR}"
    ARRAY_INDEX=$((NUM_BATCHES - 1))
    log "Submitting ${NUM_BATCHES} independent batch-chain jobs with concurrency ${CHAIN_CONCURRENT_BATCHES}"
    sbatch --array=0-${ARRAY_INDEX}%${CHAIN_CONCURRENT_BATCHES} "$0"
    exit 0
fi

mapfile -t BATCH_FILES < <(find "${BATCH_LIST_DIR}" -maxdepth 1 -type f -name 'sample_batch_*' | sort)
BATCH_FILE="${BATCH_FILES[${SLURM_ARRAY_TASK_ID}]:-}"
[[ -n "${BATCH_FILE}" ]] || fail "Could not find batch file for task ID ${SLURM_ARRAY_TASK_ID} in ${BATCH_LIST_DIR}"
run_chain "${BATCH_FILE}" "${SLURM_ARRAY_TASK_ID}"
