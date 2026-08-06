#!/usr/bin/env bash
#SBATCH --job-name=NF_Primate_Chain
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --requeue
#SBATCH --signal=B:USR1@300
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
REF_DIR="${REF_DIR:-}"
GLOBAL_REF_DIR="${GLOBAL_REF_DIR:-}"
NUCLEAR_ONLY_REF_DIR="${NUCLEAR_ONLY_REF_DIR:-}"
BATCH_SIZE="${BATCH_SIZE:-10}"
CHAIN_CONCURRENT_BATCHES="${CHAIN_CONCURRENT_BATCHES:-3}"
NUMT_CONCURRENT="${NUMT_CONCURRENT:-${CONCURRENT:-${CHAIN_CONCURRENT_BATCHES}}}"
POLL_SECONDS="${POLL_SECONDS:-120}"
LOG_DIR="${LOG_DIR:-log_all_per_batch}"
BATCH_LIST_DIR="${BATCH_LIST_DIR:-${NF_BASE_WORK_DIR}/batch_lists}"
NF_CONFIG_FILE="${NF_CONFIG_FILE:-nextflow.config}"
CLEAN_ON_SUCCESS="${CLEAN_ON_SUCCESS:-1}"
ENABLE_CHUNKED_ALIGNMENT="${ENABLE_CHUNKED_ALIGNMENT:-true}"
NEXTFLOW_MODULE="${NEXTFLOW_MODULE:-}"
export NEXTFLOW_MODULE

ACTIVE_CHILD_PID=""
ACTIVE_STAGE=""
ACTIVE_BATCH_NAME=""
ACTIVE_BATCH_FILE=""
ACTIVE_STAGE_LOG=""
REQUEUE_IN_PROGRESS=0

CLASSIFIED_FAILURE_CLASS="UNKNOWN"
CLASSIFIED_FAILURE_REASON="UNKNOWN"
CLASSIFIED_RESUME_RECOMMENDED=0
CLASSIFIED_DEFER_RECOMMENDED=1
CLASSIFIED_CLEANUP_AFTER_TERMINAL=1

log() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" >&2; }
fail() { echo "ERROR: $*" >&2; exit 1; }

current_slurm_element_id() {
    if [[ -n "${SLURM_ARRAY_JOB_ID:-}" &&
          -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        printf '%s_%s\n' \
            "${SLURM_ARRAY_JOB_ID}" \
            "${SLURM_ARRAY_TASK_ID}"
    else
        printf '%s\n' "${SLURM_JOB_ID:-}"
    fi
}

status_time() { date -u +'%Y-%m-%dT%H:%M:%SZ'; }

write_batch_status_marker() {
    local marker_kind="$1" stage="$2" status="$3" exit_code="$4"
    local failure_class="$5" failure_reason="$6" resume="$7" defer="$8" cleanup="$9"
    local stage_log="${10:-}" failed_samples="${11:-}" notes="${12:-}"
    local status_dir="${NF_BASE_WORK_DIR}/batch_status"
    local target="${status_dir}/${ACTIVE_BATCH_NAME}.${marker_kind}.tsv" tmp
    mkdir -p "${status_dir}"
    tmp="$(mktemp "${status_dir}/.${ACTIVE_BATCH_NAME}.${marker_kind}.XXXXXX")"
    # One header and one record makes markers both human-readable and machine-safe.
    printf '%b\n' 'wave_work_root\tbatch_name\tbatch_file\tstage\tstatus\texit_code\tfailure_class\tfailure_reason\tresume_recommended\tdefer_recommended\tcleanup_after_terminal\tslurm_job_id\tslurm_array_job_id\tslurm_array_task_id\tstart_time\tend_time\tstage_log\tfailed_samples_file\tnotes' > "${tmp}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${NF_BASE_WORK_DIR}" "${ACTIVE_BATCH_NAME}" "${ACTIVE_BATCH_FILE}" "$stage" "$status" "$exit_code" \
        "$failure_class" "$failure_reason" "$resume" "$defer" "$cleanup" "${SLURM_JOB_ID:-}" \
        "${SLURM_ARRAY_JOB_ID:-}" "${SLURM_ARRAY_TASK_ID:-}" "${BATCH_START_TIME:-$(status_time)}" \
        "$(status_time)" "$stage_log" "$failed_samples" "${notes//$'\t'/ }" >> "${tmp}"
    mv -f "${tmp}" "${target}"
    log "INFO: failure marker written=${target}"
}

write_batch_complete_marker() {
    write_batch_status_marker complete validation COMPLETE 0 NONE NONE 0 0 1 "${ACTIVE_STAGE_LOG}" "" "all batch outputs validated"
    rm -f "${NF_BASE_WORK_DIR}/batch_status/${ACTIVE_BATCH_NAME}.failure.tsv"
}

write_batch_failure_marker() {
    local rc="$1" failed_samples="${2:-}" notes="${3:-terminal incomplete batch may be deferred by manager}"
    if [[ -z "$failed_samples" && -s "${NF_BASE_WORK_DIR}/batch_status/${ACTIVE_BATCH_NAME}.failed_samples.tsv" ]]; then
        failed_samples="${NF_BASE_WORK_DIR}/batch_status/${ACTIVE_BATCH_NAME}.failed_samples.tsv"
    fi
    rm -f "${NF_BASE_WORK_DIR}/batch_status/${ACTIVE_BATCH_NAME}.complete.tsv"
    write_batch_status_marker failure "${ACTIVE_STAGE:-unknown}" FAILED "$rc" \
        "${CLASSIFIED_FAILURE_CLASS}" "${CLASSIFIED_FAILURE_REASON}" \
        "${CLASSIFIED_RESUME_RECOMMENDED}" "${CLASSIFIED_DEFER_RECOMMENDED}" \
        "${CLASSIFIED_CLEANUP_AFTER_TERMINAL}" "${ACTIVE_STAGE_LOG}" "$failed_samples" "$notes"
    log "INFO: terminal incomplete batch may be deferred by manager"
}

classify_stage_failure() {
    local stage="$1" rc="$2" log_file="$3" states="" text=""
    [[ -r "$log_file" ]] && text="$(cat "$log_file")"
    if command -v sacct >/dev/null 2>&1 && [[ -n "${SLURM_JOB_ID:-}" ]]; then
        states="$(sacct -j "${SLURM_JOB_ID}" --format=State --noheader 2>/dev/null || true)"
    fi
    CLASSIFIED_FAILURE_CLASS=UNKNOWN; CLASSIFIED_FAILURE_REASON=UNKNOWN
    CLASSIFIED_RESUME_RECOMMENDED=0; CLASSIFIED_DEFER_RECOMMENDED=1; CLASSIFIED_CLEANUP_AFTER_TERMINAL=1
    if [[ "$text" =~ CRAM_VALIDATION_FAILED ]]; then
        CLASSIFIED_FAILURE_CLASS=OUTPUT_INCOMPLETE; CLASSIFIED_FAILURE_REASON=CRAM_VALIDATION_FAILED
    elif [[ "$stage" == requeue || "$text" =~ (USR1|[Rr]equeue) ]]; then
        CLASSIFIED_FAILURE_CLASS=INFRASTRUCTURE; CLASSIFIED_FAILURE_REASON=TIMEOUT_SIGNAL
        CLASSIFIED_RESUME_RECOMMENDED=1; CLASSIFIED_DEFER_RECOMMENDED=0; CLASSIFIED_CLEANUP_AFTER_TERMINAL=0
    elif [[ "$states" =~ (TIMEOUT|PREEMPTED|NODE_FAIL|BOOT_FAIL) ]]; then
        CLASSIFIED_FAILURE_CLASS=INFRASTRUCTURE; CLASSIFIED_FAILURE_REASON="${BASH_REMATCH[1]}"
        CLASSIFIED_RESUME_RECOMMENDED=1; CLASSIFIED_DEFER_RECOMMENDED=0; CLASSIFIED_CLEANUP_AFTER_TERMINAL=0
    elif [[ "$text" =~ No\ space\ left\ on\ device ]]; then
        CLASSIFIED_FAILURE_CLASS=TRANSIENT_IO; CLASSIFIED_FAILURE_REASON=NO_SPACE; CLASSIFIED_RESUME_RECOMMENDED=1
    elif [[ "$text" =~ Disk\ quota\ exceeded ]]; then
        CLASSIFIED_FAILURE_CLASS=TRANSIENT_IO; CLASSIFIED_FAILURE_REASON=QUOTA_EXCEEDED; CLASSIFIED_RESUME_RECOMMENDED=1
    elif [[ "$text" =~ (Input/output\ error|Transport\ endpoint\ is\ not\ connected) ]]; then
        CLASSIFIED_FAILURE_CLASS=TRANSIENT_IO; CLASSIFIED_FAILURE_REASON=IO_ERROR; CLASSIFIED_RESUME_RECOMMENDED=1
    elif [[ "$text" =~ Read-only\ file\ system ]]; then
        CLASSIFIED_FAILURE_CLASS=TRANSIENT_IO; CLASSIFIED_FAILURE_REASON=READ_ONLY_FILESYSTEM; CLASSIFIED_RESUME_RECOMMENDED=1
    elif [[ "$text" =~ Stale\ file\ handle ]]; then
        CLASSIFIED_FAILURE_CLASS=TRANSIENT_IO; CLASSIFIED_FAILURE_REASON=STALE_FILE_HANDLE; CLASSIFIED_RESUME_RECOMMENDED=1
    elif [[ "$text" =~ (unexpected\ end\ of\ file|CRC\ error|invalid\ compressed\ data|truncated\ FASTQ|corrupt\ FASTQ|invalid\ CRAM|samtools\ quickcheck\ failed) ]]; then
        CLASSIFIED_FAILURE_CLASS=BAD_INPUT; CLASSIFIED_FAILURE_REASON=CORRUPT_INPUT
    elif [[ "$text" =~ (WGS\ ref\ missing|missing\ reference|reference\ not\ found|unsupported\ species|unsupported\ reference|malformed\ sample\ row|invalid\ metadata) ]]; then
        CLASSIFIED_FAILURE_CLASS=DETERMINISTIC; CLASSIFIED_FAILURE_REASON=INVALID_CONFIGURATION
    fi
    log "INFO: failure_class=${CLASSIFIED_FAILURE_CLASS}"
    log "INFO: failure_reason=${CLASSIFIED_FAILURE_REASON}"
    log "INFO: resume_recommended=${CLASSIFIED_RESUME_RECOMMENDED}"
    log "INFO: defer_recommended=${CLASSIFIED_DEFER_RECOMMENDED}"
    log "INFO: cleanup_after_terminal=${CLASSIFIED_CLEANUP_AFTER_TERMINAL}"
}

run_child_stage() {
    local stage_name="$1" stage_log="$2"
    shift 2

    ACTIVE_STAGE="${stage_name}"
    ACTIVE_STAGE_LOG="${stage_log}"
    mkdir -p "$(dirname "${stage_log}")"
    log "INFO: ACTIVE_STAGE=${ACTIVE_STAGE}"

    set +e
    "$@" > >(tee -a "${stage_log}") 2> >(tee -a "${stage_log}" >&2) &
    ACTIVE_CHILD_PID=$!
    wait "${ACTIVE_CHILD_PID}"
    local rc=$?
    set -e

    ACTIVE_CHILD_PID=""
    if (( rc != 0 )); then
        classify_stage_failure "${stage_name}" "${rc}" "${stage_log}"
        write_batch_failure_marker "${rc}"
    fi
    ACTIVE_STAGE=""; ACTIVE_STAGE_LOG=""
    return "${rc}"
}

handle_chain_requeue() {
    if [[ "${REQUEUE_IN_PROGRESS}" == 1 ]]; then
        return
    fi
    REQUEUE_IN_PROGRESS=1

    local element_id
    element_id="$(current_slurm_element_id)"

    log "INFO: Received USR1 before walltime"
    log "INFO: ACTIVE_BATCH_NAME=${ACTIVE_BATCH_NAME:-none}"
    log "INFO: ACTIVE_STAGE=${ACTIVE_STAGE:-none}"
    log "INFO: Preserving NF_BASE_WORK_DIR=${NF_BASE_WORK_DIR}"
    log "INFO: Requeueing array element=${element_id}"
    if [[ -n "${ACTIVE_BATCH_NAME}" ]]; then
        ACTIVE_STAGE=requeue
        CLASSIFIED_FAILURE_CLASS=INFRASTRUCTURE; CLASSIFIED_FAILURE_REASON=TIMEOUT_SIGNAL
        CLASSIFIED_RESUME_RECOMMENDED=1; CLASSIFIED_DEFER_RECOMMENDED=0; CLASSIFIED_CLEANUP_AFTER_TERMINAL=0
        rm -f "${NF_BASE_WORK_DIR}/batch_status/${ACTIVE_BATCH_NAME}.complete.tsv"
        write_batch_status_marker failure requeue REQUEUE_REQUESTED 0 INFRASTRUCTURE TIMEOUT_SIGNAL 1 0 0 "${ACTIVE_STAGE_LOG}" "" "self-requeue requested"
        log "INFO: preserving work during self-requeue"
    fi

    if [[ -n "${ACTIVE_CHILD_PID}" ]] &&
       kill -0 "${ACTIVE_CHILD_PID}" 2>/dev/null; then

        kill -TERM "${ACTIVE_CHILD_PID}" 2>/dev/null || true

        for _ in $(seq 1 60); do
            kill -0 "${ACTIVE_CHILD_PID}" 2>/dev/null || break
            sleep 2
        done

        if kill -0 "${ACTIVE_CHILD_PID}" 2>/dev/null; then
            kill -KILL "${ACTIVE_CHILD_PID}" 2>/dev/null || true
        fi

        wait "${ACTIVE_CHILD_PID}" 2>/dev/null || true
        ACTIVE_CHILD_PID=""
    fi

    [[ -n "${element_id}" ]] || {
        log "ERROR: Cannot determine Slurm job/array element ID"
        exit 75
    }

    scontrol requeue "${element_id}" || {
        log "ERROR: Failed to requeue ${element_id}"
        exit 75
    }

    exit 0
}

trap handle_chain_requeue USR1

load_nextflow() {
    if command -v nextflow >/dev/null 2>&1; then
        echo "INFO: Using existing Nextflow: $(command -v nextflow)" >&2
        nextflow -version
        return 0
    fi

    if [[ -z "${NEXTFLOW_MODULE:-}" ]]; then
        echo "ERROR: NEXTFLOW_MODULE is not set and nextflow is not already available" >&2
        module spider Nextflow >&2 || true
        exit 1
    fi

    echo "INFO: Loading Nextflow module: ${NEXTFLOW_MODULE}" >&2
    module load "${NEXTFLOW_MODULE}" || {
        echo "ERROR: Failed to load ${NEXTFLOW_MODULE}" >&2
        module spider Nextflow >&2 || true
        exit 1
    }

    command -v nextflow >/dev/null 2>&1 || {
        echo "ERROR: nextflow not found after loading ${NEXTFLOW_MODULE}" >&2
        exit 1
    }

    echo "INFO: Nextflow executable: $(command -v nextflow)" >&2
    nextflow -version
}

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
echo "INFO: ENABLE_CHUNKED_ALIGNMENT=${ENABLE_CHUNKED_ALIGNMENT}"
echo "INFO: NEXTFLOW_MODULE=${NEXTFLOW_MODULE:-<none>}"

cleanup_work_dir_if_requested() {
    local stage_name="$1"
    local work_dir="$2"

    if [[ "${CLEAN_ON_SUCCESS}" == "1" || "${CLEAN_ON_SUCCESS}" == "true" ]]; then
        local base_real target_real rel
        base_real="$(readlink -f "${NF_BASE_WORK_DIR}")"
        target_real="$(readlink -f -m "${work_dir}")"
        rel="${target_real#"${base_real}"/}"
        if [[ "${target_real}" == "/" || "${target_real}" == "${base_real}" || "${target_real}" != "${base_real}/"* ||
              ! "${rel}" =~ ^(pre|numt|round1|round2)/batch_[A-Za-z0-9._-]+$ ]]; then
            log "ERROR: refusing unsafe batch cache cleanup: ${work_dir} -> ${target_real}"
            return 64
        fi
        if [[ -d "${target_real}" ]]; then
            log "CLEAN_ON_SUCCESS=${CLEAN_ON_SUCCESS}; deleting ${stage_name} work directory: ${work_dir}"
            rm -rf -- "${target_real}"
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

    load_nextflow

    local numt_cmd=(
        nextflow
        -C "${NF_CONFIG_PATH}"
        run "${NUMT_LAUNCH_SCRIPT}"
        -profile cluster
        -work-dir "${numt_batch_work_dir}"
        -with-report "${numt_batch_dir}/nextflow.report.html"
        -with-timeline "${numt_batch_dir}/nextflow.timeline.html"
        -with-trace "${numt_batch_dir}/nextflow.trace.tsv"
        -queue-size "${NUMT_CONCURRENT}"
        --numt_config "${auto_config}"
    )

    run_child_stage numt "${LOG_DIR}/${batch_name}.numt.log" env \
        NF_CONFIG_PATH="${NF_CONFIG_PATH}" \
        NXF_OPTS="${NXF_OPTS:-}" \
        "${numt_cmd[@]}"
}

validate_pre_to_round1() {
    local batch_file="$1" missing=0 sample cram crai
    log "Validating preprocessing outputs for ${batch_file}"
    while IFS= read -r sample; do
        cram="${PRE_OUTPUT_DIR}/${sample}/alignment/${sample}.cram"; crai="${cram}.crai"
        [[ -s "${cram}" ]] || { echo "MISSING/EMPTY: ${cram}" >&2; missing=1; }
        [[ -s "${crai}" ]] || { echo "MISSING/EMPTY: ${crai}" >&2; missing=1; }
    done < <(read_sample_ids "${batch_file}")
    [[ "${missing}" -eq 0 ]]
}

write_failed_samples() {
    local batch_file="$1" output="${NF_BASE_WORK_DIR}/batch_status/${ACTIVE_BATCH_NAME}.failed_samples.tsv"
    local tmp sample cram crai now
    mkdir -p "$(dirname "$output")"
    tmp="$(mktemp "$(dirname "$output")/.${ACTIVE_BATCH_NAME}.failed_samples.XXXXXX")"
    printf 'sample_id\tfailure_reason\texpected_cram\texpected_crai\tvalidation_time\n' > "$tmp"
    now="$(status_time)"
    while IFS= read -r sample; do
        cram="${PRE_OUTPUT_DIR}/${sample}/alignment/${sample}.cram"; crai="${cram}.crai"
        if [[ ! -s "$cram" || ! -s "$crai" ]]; then
            printf '%s\t%s\t%s\t%s\t%s\n' "$sample" CRAM_VALIDATION_FAILED "$cram" "$crai" "$now" >> "$tmp"
        fi
    done < <(read_sample_ids "$batch_file")
    mv -f "$tmp" "$output"
    log "INFO: failed samples file=${output}"
}

record_incomplete_outputs() {
    local reason="$1" rc="${2:-1}"
    ACTIVE_STAGE=validation
    CLASSIFIED_FAILURE_CLASS=OUTPUT_INCOMPLETE; CLASSIFIED_FAILURE_REASON="$reason"
    CLASSIFIED_RESUME_RECOMMENDED=0; CLASSIFIED_DEFER_RECOMMENDED=1; CLASSIFIED_CLEANUP_AFTER_TERMINAL=1
    write_batch_failure_marker "$rc"
    ACTIVE_STAGE=""
    return "$rc"
}

validate_numt_to_round1() {
    local batch_file="$1" missing=0 sample bed
    log "Validating NUMT outputs for ${batch_file}"
    while IFS= read -r sample; do
        bed="${NUMT_BESTHIT_OUTDIR}/${sample}.highconf_numt.bed"
        [[ -e "${bed}" ]] || { echo "MISSING: ${bed}" >&2; missing=1; }
    done < <(read_sample_ids "${batch_file}")
    [[ "${missing}" -eq 0 ]]
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
    [[ "${missing}" -eq 0 ]]
}

validate_round2_final() {
    local batch_file="$1" missing=0 sample out_dir
    log "Validating round 2 outputs for ${batch_file}"
    while IFS= read -r sample; do
        out_dir="${ROUND_OUTPUT_DIR}/${sample}/round_2_variant_calling_original_coords"
        [[ -d "${out_dir}" ]] || { echo "MISSING DIR: ${out_dir}" >&2; missing=1; continue; }
        find "${out_dir}" -type f -name "${sample}.round2.original_coords.clean.final.split.vcf.gz" -size +0c -print -quit 2>/dev/null | grep -q . || { echo "MISSING/EMPTY round2 final VCF under: ${out_dir}" >&2; missing=1; }
    done < <(read_sample_ids "${batch_file}")
    [[ "${missing}" -eq 0 ]]
}

run_chain() {
    local batch_file="$1" batch_id="$2"
    local batch_name="batch_${batch_id}"
    local pre_batch_work_dir="${PRE_NF_BASE_WORK_DIR}/${batch_name}"
    local round1_batch_work_dir="${ROUND1_NF_BASE_WORK_DIR}/${batch_name}"
    local round2_batch_work_dir="${ROUND2_NF_BASE_WORK_DIR}/${batch_name}"

    ACTIVE_BATCH_NAME="${batch_name}"
    ACTIVE_BATCH_FILE="$(readlink -f "${batch_file}")"
    BATCH_START_TIME="$(status_time)"
    log "INFO: ACTIVE_BATCH_NAME=${ACTIVE_BATCH_NAME}"
    write_batch_status_marker failure initialization RUNNING 0 NONE NONE 1 0 0 "" "" "batch chain started"

    log "Starting per-batch chain ${batch_id}: ${batch_file}"
    run_child_stage pre "${LOG_DIR}/${batch_name}.pre.log" env \
        ENABLE_CHUNKED_ALIGNMENT="${ENABLE_CHUNKED_ALIGNMENT}" \
        BATCH_FILE="${batch_file}" \
        BATCH_ID="${batch_name}" \
        FULL_SAMPLE_LIST="${batch_file}" \
        OUTPUT_DIR="${PRE_OUTPUT_DIR}" \
        NF_BASE_WORK_DIR="${PRE_NF_BASE_WORK_DIR}" \
        NF_CONFIG_FILE="${NF_CONFIG_FILE}" \
        NEXTFLOW_MODULE="${NEXTFLOW_MODULE}" \
        DEFER_WORK_DIR_CLEANUP=1 \
        CHAIN_MANAGED_REQUEUE=1 \
        FAILED_SAMPLES_FILE="${NF_BASE_WORK_DIR}/batch_status/${batch_name}.failed_samples.tsv" \
        bash "${PRE_LAUNCH_SCRIPT}"
    if ! validate_pre_to_round1 "${batch_file}"; then
        write_failed_samples "${batch_file}"
        ACTIVE_STAGE=validation; ACTIVE_STAGE_LOG="${LOG_DIR}/${batch_name}.pre.log"
        CLASSIFIED_FAILURE_CLASS=OUTPUT_INCOMPLETE; CLASSIFIED_FAILURE_REASON=CRAM_VALIDATION_FAILED
        CLASSIFIED_RESUME_RECOMMENDED=0; CLASSIFIED_DEFER_RECOMMENDED=1; CLASSIFIED_CLEANUP_AFTER_TERMINAL=1
        write_batch_failure_marker 1 "${NF_BASE_WORK_DIR}/batch_status/${batch_name}.failed_samples.tsv"
        return 1
    fi

    if validate_numt_to_round1 "${batch_file}"; then
        log "Skipping NUMT for ${batch_id}; expected NUMT outputs already exist."
    else
        log "NUMT outputs incomplete for ${batch_id}; launching NUMT step."
        run_numt_nextflow "${batch_file}" "${batch_id}"
        if ! validate_numt_to_round1 "${batch_file}"; then
            record_incomplete_outputs NUMT_VALIDATION_FAILED
            return 1
        fi
    fi

    if validate_round1_to_round2 "${batch_file}"; then
        log "Skipping round1 for ${batch_id}; expected round1 outputs already exist."
    else
        log "Round1 outputs incomplete for ${batch_id}; launching round1 step."
        log "Launching round1 for ${batch_id}"
        log "  BATCH_FILE=${batch_file}"
        log "  OUTPUT_DIR=${ROUND1_OUTDIR}"
        log "  CRAM_DIRS=${PRE_OUTPUT_DIR}"
        log "  NUMT_BED_DIR=${NUMT_BESTHIT_OUTDIR}"
        log "  NF_BASE_WORK_DIR=${ROUND1_NF_BASE_WORK_DIR}"
        run_child_stage round1 "${LOG_DIR}/${batch_name}.round1.log" env \
            BATCH_FILE="${batch_file}" \
            BATCH_ID="${batch_name}" \
            FULL_SAMPLE_LIST="${batch_file}" \
            OUTPUT_DIR="${ROUND1_OUTDIR}" \
            CRAM_DIRS="${PRE_OUTPUT_DIR}" \
            NUMT_BED_DIR="${NUMT_BESTHIT_OUTDIR}" \
            NF_BASE_WORK_DIR="${ROUND1_NF_BASE_WORK_DIR}" \
            NF_CONFIG_FILE="${NF_CONFIG_FILE}" \
            NEXTFLOW_MODULE="${NEXTFLOW_MODULE}" \
            DEFER_WORK_DIR_CLEANUP=1 \
            CHAIN_MANAGED_REQUEUE=1 \
            bash "${ROUND1_LAUNCH_SCRIPT}"
        if ! validate_round1_to_round2 "${batch_file}"; then
            record_incomplete_outputs ROUND1_VALIDATION_FAILED
            return 1
        fi
    fi

    if validate_round2_final "${batch_file}"; then
        log "Skipping round2 for ${batch_id}; expected round2 outputs already exist."
    else
        log "Round2 outputs incomplete for ${batch_id}; launching round2 step."
        log "Launching round2 for ${batch_id}"
        log "  BATCH_FILE=${batch_file}"
        log "  OUTPUT_DIR=${ROUND_OUTPUT_DIR}"
        log "  ROUND1_OUTDIR=${ROUND1_OUTDIR}"
        log "  NF_BASE_WORK_DIR=${ROUND2_NF_BASE_WORK_DIR}"
        run_child_stage round2 "${LOG_DIR}/${batch_name}.round2.log" env \
            BATCH_FILE="${batch_file}" \
            BATCH_ID="${batch_name}" \
            FULL_SAMPLE_LIST="${batch_file}" \
            OUTPUT_DIR="${ROUND_OUTPUT_DIR}" \
            ROUND1_OUTDIR="${ROUND1_OUTDIR}" \
            NF_BASE_WORK_DIR="${ROUND2_NF_BASE_WORK_DIR}" \
            NF_CONFIG_FILE="${NF_CONFIG_FILE}" \
            NEXTFLOW_MODULE="${NEXTFLOW_MODULE}" \
            DEFER_WORK_DIR_CLEANUP=1 \
            CHAIN_MANAGED_REQUEUE=1 \
            bash "${ROUND2_LAUNCH_SCRIPT}"
        if ! validate_round2_final "${batch_file}"; then
            record_incomplete_outputs ROUND2_VALIDATION_FAILED
            return 1
        fi
    fi
    write_batch_complete_marker
    cleanup_work_dir_if_requested "pre_${batch_name}" "${pre_batch_work_dir}"
    cleanup_work_dir_if_requested "numt_${batch_name}" "${NUMT_NF_BASE_WORK_DIR}/${batch_name}"
    cleanup_work_dir_if_requested "round1_${batch_name}" "${round1_batch_work_dir}"
    cleanup_work_dir_if_requested "round2_${batch_name}" "${round2_batch_work_dir}"
    log "INFO: successful batch cache cleanup completed"
    log "Completed per-batch chain ${batch_id}"
    ACTIVE_BATCH_NAME=""; ACTIVE_BATCH_FILE=""; ACTIVE_STAGE=""; ACTIVE_STAGE_LOG=""
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
