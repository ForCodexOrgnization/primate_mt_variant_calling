#!/usr/bin/env bash
#SBATCH --job-name=NF_Primate_All
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=log_all/nf_all_%j.log

set -euo pipefail

# ==============================================================================
# Chain preprocessing -> NUMT discovery -> round 1 -> round 2 and validate hand-off files.
# Run from the repository directory on a login/submission node:
#   bash launch_pipeline_all.sh
# or submit the coordinator itself:
#   sbatch launch_pipeline_all.sh
#
# All values can be overridden from the environment, for example:
#   FULL_SAMPLE_LIST=/path/samples.tsv PRE_OUTPUT_DIR=/path/pre ROUND_OUTPUT_DIR=/path/results bash launch_pipeline_all.sh
# ==============================================================================

FULL_SAMPLE_LIST="${FULL_SAMPLE_LIST:-/nfs/roberts/project/pi_njl27/lt692/primate_mt_variant_calling/human_sample.txt}"
PRE_OUTPUT_DIR="${PRE_OUTPUT_DIR:-/nfs/roberts/scratch/pi_njl27/lt692/primate_results_test}"
ROUND_OUTPUT_DIR="${ROUND_OUTPUT_DIR:-/nfs/roberts/scratch/pi_njl27/lt692/primate_results_test}"
ROUND1_OUTDIR="${ROUND1_OUTDIR:-${ROUND_OUTPUT_DIR}}"

PRE_LAUNCH_SCRIPT="${PRE_LAUNCH_SCRIPT:-launch_pipeline_pre.sh}"
NUMT_LAUNCH_SCRIPT="${NUMT_LAUNCH_SCRIPT:-launch_pipeline_numt.sh}"
ROUND1_LAUNCH_SCRIPT="${ROUND1_LAUNCH_SCRIPT:-launch_pipeline_round1.sh}"
ROUND2_LAUNCH_SCRIPT="${ROUND2_LAUNCH_SCRIPT:-launch_pipeline_round2.sh}"

NUMT_DISCOVERY_OUTROOT="${NUMT_DISCOVERY_OUTROOT:-${ROUND_OUTPUT_DIR}/numt_discovery}"
NUMT_BESTHIT_OUTDIR="${NUMT_BESTHIT_OUTDIR:-${ROUND_OUTPUT_DIR}/numt_besthit}"
# NUMT discovery reference inputs. Use the same naming as nextflow.config:
#   GLOBAL_REF_DIR == whole-genome reference directory
#   REF_DIR        == chrM reference directory
# NUCLEAR_ONLY_REF_DIR is only used by numt_detection and must point to a
# chrM-excluded reference set. Round 1 only consumes the resulting
# high-confidence BEDs via NUMT_BED_DIR.
REF_DIR="${REF_DIR:-/home/lt692/scratch_pi_njl27/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_chrM}"
GLOBAL_REF_DIR="${GLOBAL_REF_DIR:-/home/lt692/scratch_pi_njl27/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_whole}"
NUCLEAR_ONLY_REF_DIR="${NUCLEAR_ONLY_REF_DIR:-/home/lt692/scratch_pi_njl27/lt692/primate_mtDNA_analysis/references/variant_calling/Ref_nuclear_only}"
NUMT_CONCURRENT="${NUMT_CONCURRENT:-${CONCURRENT:-2}}"

POLL_SECONDS="${POLL_SECONDS:-120}"
LOG_DIR="${LOG_DIR:-log_all}"
mkdir -p "${LOG_DIR}"

# When this coordinator is submitted with sbatch, Slurm may execute a spool copy of
# the script. Prefer the original submission directory so relative launch script paths
# still resolve to the repository checkout that contains launch_pipeline_pre.sh, etc.
if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/preprocessing.nf" ]]; then
    repo_dir="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
else
    repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
cd "${repo_dir}"

log() {
    printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" >&2
}

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

read_sample_ids() {
    awk -F '\t' 'NF >= 2 && $1 !~ /^#/ && $1 != "" {print $1}' "${FULL_SAMPLE_LIST}"
}

submit_step() {
    local step_name="$1"
    local script_path="$2"
    shift 2

    [[ -f "${script_path}" ]] || fail "Missing launch script: ${script_path}"

    local submit_log="${LOG_DIR}/${step_name}.submit.log"
    log "Submitting ${step_name} via ${script_path}"

    env "$@" bash "${script_path}" 2>&1 | tee "${submit_log}" >&2

    local job_id
    job_id="$({ sed -n 's/.*Submitted batch job \([0-9][0-9]*\).*/\1/p' "${submit_log}" || true; } | tail -n 1)"
    [[ -n "${job_id}" ]] || fail "Could not parse Slurm job id for ${step_name}; see ${submit_log}"

    log "${step_name} Slurm job id: ${job_id}"
    echo "${job_id}"
}

job_state_sacct() {
    local job_id="$1"
    sacct -j "${job_id}" --format=State --noheader 2>/dev/null \
        | awk 'NF {print $1}' \
        | sed 's/+$//' \
        | sort -u \
        | paste -sd, -
}

wait_for_job() {
    local job_id="$1"
    local step_name="$2"

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
        else
            if ! squeue -h -j "${job_id}" >/dev/null 2>&1 || [[ -z "$(squeue -h -j "${job_id}" 2>/dev/null)" ]]; then
                log "${step_name} job ${job_id} is no longer in squeue; sacct unavailable, assuming completion"
                return 0
            fi
            log "${step_name} job ${job_id} still present in squeue; sleeping ${POLL_SECONDS}s"
        fi
        sleep "${POLL_SECONDS}"
    done
}

validate_pre_to_round1() {
    log "Validating preprocessing outputs needed by round 1"
    local missing=0 sample cram crai
    while IFS= read -r sample; do
        cram="${PRE_OUTPUT_DIR}/${sample}/alignment/${sample}.cram"
        crai="${PRE_OUTPUT_DIR}/${sample}/alignment/${sample}.cram.crai"
        [[ -s "${cram}" ]] || { echo "MISSING/EMPTY: ${cram}" >&2; missing=1; }
        [[ -s "${crai}" ]] || { echo "MISSING/EMPTY: ${crai}" >&2; missing=1; }
    done < <(read_sample_ids)
    [[ "${missing}" -eq 0 ]] || fail "Preprocessing outputs are incomplete; round 1 was not submitted"
    log "Preprocessing validation passed"
}

validate_numt_to_round1() {
    log "Validating NUMT discovery outputs needed by round 1"
    local missing=0 sample bed
    while IFS= read -r sample; do
        bed="${NUMT_BESTHIT_OUTDIR}/${sample}.highconf_numt.bed"
        # Empty BED files are valid for samples with no high-confidence NUMT calls.
        [[ -e "${bed}" ]] || { echo "MISSING: ${bed}" >&2; missing=1; }
    done < <(read_sample_ids)
    [[ "${missing}" -eq 0 ]] || fail "NUMT discovery outputs are incomplete; round 1 was not submitted"
    log "NUMT discovery validation passed"
}

numt_outputs_complete() {
    local missing=0 sample bed
    while IFS= read -r sample; do
        bed="${NUMT_BESTHIT_OUTDIR}/${sample}.highconf_numt.bed"
        [[ -e "${bed}" ]] || missing=1
    done < <(read_sample_ids)
    [[ "${missing}" -eq 0 ]]
}

validate_round1_to_round2() {
    log "Validating round 1 outputs needed by round 2"
    local missing=0 sample
    while IFS= read -r sample; do
        local bam="${ROUND1_OUTDIR}/${sample}/round_1/candidate_reads/${sample}.with_mates.bam"
        local vcf_dir="${ROUND1_OUTDIR}/${sample}/round_1_variant_calling_decoy"
        local numt_fa="${ROUND1_OUTDIR}/${sample}/round_1/numt_decoy_ref/${sample}.original_numt.fa"
        local numt_vcf="${ROUND1_OUTDIR}/${sample}/round_1/numt_decoy_variant_calling/${sample}.numt_decoy.raw.vcf.gz"
        local numt_tbi="${numt_vcf}.tbi"
        [[ -s "${bam}" ]] || { echo "MISSING/EMPTY: ${bam}" >&2; missing=1; }
        [[ -d "${vcf_dir}" ]] || { echo "MISSING DIR: ${vcf_dir}" >&2; missing=1; }
        find "${vcf_dir}" -type f -name "${sample}.numt_decoy.clean.final.split.vcf" -size +0c -print -quit 2>/dev/null | grep -q . || { echo "MISSING/EMPTY round1 final VCF under: ${vcf_dir}" >&2; missing=1; }
        [[ -e "${numt_fa}" ]] || { echo "MISSING: ${numt_fa}" >&2; missing=1; }
        [[ -s "${numt_vcf}" ]] || { echo "MISSING/EMPTY: ${numt_vcf}" >&2; missing=1; }
        [[ -s "${numt_tbi}" ]] || { echo "MISSING/EMPTY: ${numt_tbi}" >&2; missing=1; }
    done < <(read_sample_ids)
    [[ "${missing}" -eq 0 ]] || fail "Round 1 outputs are incomplete; round 2 was not submitted"
    log "Round 1 validation passed"
}

validate_round2_final() {
    log "Validating round 2 final outputs"
    local missing=0 sample out_dir
    while IFS= read -r sample; do
        out_dir="${ROUND_OUTPUT_DIR}/${sample}/round_2_variant_calling_original_coords"
        [[ -d "${out_dir}" ]] || { echo "MISSING DIR: ${out_dir}" >&2; missing=1; continue; }
        find "${out_dir}" -type f -name "${sample}.round2.original_coords.clean.final.split.vcf.gz" -size +0c -print -quit 2>/dev/null | grep -q . || { echo "MISSING/EMPTY round2 final VCF under: ${out_dir}" >&2; missing=1; }
    done < <(read_sample_ids)
    [[ "${missing}" -eq 0 ]] || fail "Round 2 outputs are incomplete"
    log "Round 2 validation passed"
}

log "Starting chained pipeline"
log "Sample list: ${FULL_SAMPLE_LIST}"
log "Preprocessing output dir: ${PRE_OUTPUT_DIR}"
log "Round output dir: ${ROUND_OUTPUT_DIR}"
log "NUMT discovery outroot: ${NUMT_DISCOVERY_OUTROOT}"
log "NUMT best-hit BED dir: ${NUMT_BESTHIT_OUTDIR}"
log "NUMT global reference dir: ${GLOBAL_REF_DIR}"
log "NUMT nuclear-only reference dir: ${NUCLEAR_ONLY_REF_DIR}"
log "NUMT chrM reference dir: ${REF_DIR}"

pre_job="$(submit_step preprocess "${PRE_LAUNCH_SCRIPT}" FULL_SAMPLE_LIST="${FULL_SAMPLE_LIST}" OUTPUT_DIR="${PRE_OUTPUT_DIR}")"
wait_for_job "${pre_job}" preprocess
validate_pre_to_round1

if numt_outputs_complete; then
    log "Skipping NUMT discovery; all high-confidence NUMT BED outputs already exist under ${NUMT_BESTHIT_OUTDIR}"
else
    numt_job="$(submit_step numt "${NUMT_LAUNCH_SCRIPT}" FULL_SAMPLE_LIST="${FULL_SAMPLE_LIST}" PRE_OUTPUT_DIR="${PRE_OUTPUT_DIR}" GLOBAL_REF_DIR="${GLOBAL_REF_DIR}" NUCLEAR_ONLY_REF_DIR="${NUCLEAR_ONLY_REF_DIR}" REF_DIR="${REF_DIR}" DISCOVERY_OUTROOT="${NUMT_DISCOVERY_OUTROOT}" BESTHIT_OUTDIR="${NUMT_BESTHIT_OUTDIR}" CONCURRENT="${NUMT_CONCURRENT}")"
    wait_for_job "${numt_job}" numt
fi
validate_numt_to_round1

round1_job="$(submit_step round1 "${ROUND1_LAUNCH_SCRIPT}" FULL_SAMPLE_LIST="${FULL_SAMPLE_LIST}" OUTPUT_DIR="${ROUND1_OUTDIR}" CRAM_DIRS="${PRE_OUTPUT_DIR}" NUMT_BED_DIR="${NUMT_BESTHIT_OUTDIR}")"
wait_for_job "${round1_job}" round1
validate_round1_to_round2

round2_job="$(submit_step round2 "${ROUND2_LAUNCH_SCRIPT}" FULL_SAMPLE_LIST="${FULL_SAMPLE_LIST}" OUTPUT_DIR="${ROUND_OUTPUT_DIR}" ROUND1_OUTDIR="${ROUND1_OUTDIR}")"
wait_for_job "${round2_job}" round2
validate_round2_final

log "All four steps completed and validation checks passed"
