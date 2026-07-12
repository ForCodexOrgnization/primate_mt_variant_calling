#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Watch completed preprocessing CRAM/CRAI outputs and launch downstream pipeline
# for newly completed samples.
#
# This script scans PRE_OUTPUT_DIR for:
#   SAMPLE/alignment/SAMPLE.cram
#   SAMPLE/alignment/SAMPLE.cram.crai
#
# If both files exist, are non-empty, and samtools quickcheck passes, the script
# extracts that sample's line from MASTER_SAMPLE_LIST, writes a one-sample TSV,
# and submits LAUNCH_SCRIPT for the downstream full pipeline.
#
# Already submitted samples are tracked in:
#   ${DISPATCH_DIR}/state/submitted_samples.txt
# ==============================================================================

PRE_OUTPUT_DIR="${PRE_OUTPUT_DIR:-/home/lt692/scratch_pi_njl27/lt692/primate_results_test}"
ROUND_OUTPUT_DIR="${ROUND_OUTPUT_DIR:-${PRE_OUTPUT_DIR}}"

MASTER_SAMPLE_LIST="${MASTER_SAMPLE_LIST:-/home/lt692/project_pi_njl27/lt692/primate_mt_variant_calling/human_sample.txt}"

REPO_DIR="${REPO_DIR:-/home/lt692/project_pi_njl27/lt692/primate_mt_variant_calling}"
LAUNCH_SCRIPT="${LAUNCH_SCRIPT:-${REPO_DIR}/launch_pipeline_streaming_per_sample.sh}"

DISPATCH_DIR="${DISPATCH_DIR:-${REPO_DIR}/auto_dispatch_completed_pre}"
STATE_DIR="${DISPATCH_DIR}/state"
INPUT_DIR="${DISPATCH_DIR}/inputs"
LOG_DIR="${DISPATCH_DIR}/logs"

SUBMITTED_FILE="${STATE_DIR}/submitted_samples.txt"
LOCK_DIR="${STATE_DIR}/dispatch.lock"

WATCH_INTERVAL_SECONDS="${WATCH_INTERVAL_SECONDS:-300}"
RUN_ONCE="${RUN_ONCE:-0}"

# Downstream launcher settings
NF_BASE_WORK_DIR="${NF_BASE_WORK_DIR:-/home/lt692/ycga_work/nf_work_dir_streaming_per_sample}"
CLEAN_ON_SUCCESS="${CLEAN_ON_SUCCESS:-1}"
MAX_CONCURRENT="${MAX_CONCURRENT:-1}"
NF_CONFIG_FILE="${NF_CONFIG_FILE:-nextflow.config}"

mkdir -p "${STATE_DIR}" "${INPUT_DIR}" "${LOG_DIR}"
touch "${SUBMITTED_FILE}"

log() {
    printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" >&2
}

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

acquire_lock() {
    if mkdir "${LOCK_DIR}" 2>/dev/null; then
        trap 'rm -rf "${LOCK_DIR}"' EXIT
    else
        fail "Another dispatcher appears to be running: ${LOCK_DIR}"
    fi
}

sample_already_submitted() {
    local sample_id="$1"
    grep -Fxq "${sample_id}" "${SUBMITTED_FILE}"
}

mark_sample_submitted() {
    local sample_id="$1"
    echo "${sample_id}" >> "${SUBMITTED_FILE}"
    sort -u "${SUBMITTED_FILE}" -o "${SUBMITTED_FILE}"
}

cram_is_complete() {
    local sample_id="$1"
    local cram="${PRE_OUTPUT_DIR}/${sample_id}/alignment/${sample_id}.cram"
    local crai="${PRE_OUTPUT_DIR}/${sample_id}/alignment/${sample_id}.cram.crai"

    [[ -s "${cram}" ]] || return 1
    [[ -s "${crai}" ]] || return 1

    if command -v samtools >/dev/null 2>&1; then
        samtools quickcheck "${cram}" >/dev/null 2>&1 || return 1
    else
        log "WARN: samtools not found; using non-empty CRAM/CRAI check only for ${sample_id}"
    fi

    return 0
}

write_one_sample_tsv() {
    local sample_id="$1"
    local out_tsv="$2"

    # Match first column exactly. This supports both 2-column and 3-column sample lists.
    awk -F '\t' -v s="${sample_id}" 'BEGIN{OFS=FS} $1 == s {print; found=1} END{exit found ? 0 : 1}' \
        "${MASTER_SAMPLE_LIST}" > "${out_tsv}"
}

submit_downstream_for_sample() {
    local sample_id="$1"
    local sample_tsv="${INPUT_DIR}/${sample_id}.sample.tsv"
    local submit_log="${LOG_DIR}/${sample_id}.submit.log"

    write_one_sample_tsv "${sample_id}" "${sample_tsv}" || {
        log "WARN: ${sample_id} has complete CRAM/CRAI but was not found in MASTER_SAMPLE_LIST: ${MASTER_SAMPLE_LIST}"
        rm -f "${sample_tsv}"
        return 1
    }

    [[ -s "${sample_tsv}" ]] || {
        log "WARN: Generated empty sample TSV for ${sample_id}; skipping"
        rm -f "${sample_tsv}"
        return 1
    }

    log "Submitting downstream pipeline for ${sample_id}"
    log "One-sample TSV: ${sample_tsv}"

    # launch_pipeline_streaming_per_sample.sh expects FULL_SAMPLE_LIST and submits/uses
    # a one-sample array. MAX_CONCURRENT=1 is safest for one-sample TSV.
    sbatch \
        --export=ALL,FULL_SAMPLE_LIST="${sample_tsv}",PRE_OUTPUT_DIR="${PRE_OUTPUT_DIR}",ROUND_OUTPUT_DIR="${ROUND_OUTPUT_DIR}",NF_BASE_WORK_DIR="${NF_BASE_WORK_DIR}",CLEAN_ON_SUCCESS="${CLEAN_ON_SUCCESS}",MAX_CONCURRENT="${MAX_CONCURRENT}",NF_CONFIG_FILE="${NF_CONFIG_FILE}" \
        "${LAUNCH_SCRIPT}" > "${submit_log}"

    cat "${submit_log}" >&2
    mark_sample_submitted "${sample_id}"
}

scan_once() {
    log "Scanning PRE_OUTPUT_DIR: ${PRE_OUTPUT_DIR}"

    local n_ready=0
    local n_submitted=0
    local n_skipped=0

    shopt -s nullglob

    for sample_dir in "${PRE_OUTPUT_DIR}"/*; do
        [[ -d "${sample_dir}" ]] || continue

        local sample_id
        sample_id="$(basename "${sample_dir}")"

        # Skip non-sample helper directories if present.
        case "${sample_id}" in
            numt_discovery|numt_besthit|logs|tmp|work)
                continue
                ;;
        esac

        if sample_already_submitted "${sample_id}"; then
            ((n_skipped++)) || true
            continue
        fi

        if cram_is_complete "${sample_id}"; then
            ((n_ready++)) || true
            if submit_downstream_for_sample "${sample_id}"; then
                ((n_submitted++)) || true
            fi
        fi
    done

    log "Scan complete: ready=${n_ready}, submitted=${n_submitted}, already_submitted_or_skipped=${n_skipped}"
}

main() {
    [[ -d "${PRE_OUTPUT_DIR}" ]] || fail "PRE_OUTPUT_DIR does not exist: ${PRE_OUTPUT_DIR}"
    [[ -s "${MASTER_SAMPLE_LIST}" ]] || fail "MASTER_SAMPLE_LIST missing or empty: ${MASTER_SAMPLE_LIST}"
    [[ -f "${LAUNCH_SCRIPT}" ]] || fail "LAUNCH_SCRIPT does not exist: ${LAUNCH_SCRIPT}"

    acquire_lock

    if command -v module >/dev/null 2>&1; then
        module load SAMtools/1.21-GCC-12.2.0 >/dev/null 2>&1 || true
    fi

    while true; do
        scan_once

        if [[ "${RUN_ONCE}" == "1" || "${RUN_ONCE}" == "true" ]]; then
            break
        fi

        log "Sleeping ${WATCH_INTERVAL_SECONDS}s before next scan"
        sleep "${WATCH_INTERVAL_SECONDS}"
    done
}

main "$@"
