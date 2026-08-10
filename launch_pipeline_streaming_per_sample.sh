#!/usr/bin/env bash
#SBATCH --job-name=NF_Primate_Stream
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --requeue
#SBATCH --signal=B:USR1@300
#SBATCH --output=log_streaming/nf_streaming_%A_%a.log

set -euo pipefail

# A stable, sample-centric launcher.  Work/cache identity is SAMPLE_ID, not a
# submission "wave"; published workflow outputs remain in the existing trees.
FULL_SAMPLE_LIST="${FULL_SAMPLE_LIST:-}"
PRE_OUTPUT_DIR="${PRE_OUTPUT_DIR:-}"
ROUND_OUTPUT_DIR="${ROUND_OUTPUT_DIR:-}"
ROUND1_OUTDIR="${ROUND1_OUTDIR:-${ROUND_OUTPUT_DIR}}"
NF_BASE_WORK_DIR="${NF_BASE_WORK_DIR:-}"
NF_CONFIG_FILE="${NF_CONFIG_FILE:-}"
GLOBAL_REF_DIR="${GLOBAL_REF_DIR:-}"
REF_DIR="${REF_DIR:-}"
NUCLEAR_ONLY_REF_DIR="${NUCLEAR_ONLY_REF_DIR:-}"
NUMT_DISCOVERY_OUTROOT="${NUMT_DISCOVERY_OUTROOT:-${ROUND_OUTPUT_DIR:+${ROUND_OUTPUT_DIR}/numt_discovery}}"
NUMT_BESTHIT_OUTDIR="${NUMT_BESTHIT_OUTDIR:-${ROUND_OUTPUT_DIR:+${ROUND_OUTPUT_DIR}/numt_besthit}}"
MAX_CONCURRENT="${MAX_CONCURRENT:-10}"
STREAM_PARTITION="${STREAM_PARTITION:-}"
PIPELINE_PROFILE="${PIPELINE_PROFILE:-cluster}"
CLEAN_VALIDATED_STAGE_WORK="${CLEAN_VALIDATED_STAGE_WORK:-1}"
REMOVE_SAMPLE_ROOT_ON_SUCCESS="${REMOVE_SAMPLE_ROOT_ON_SUCCESS:-1}"
IMMEDIATE_SAMPLE_RETRIES="${IMMEDIATE_SAMPLE_RETRIES:-1}"
IMMEDIATE_RETRY_DELAY_SECONDS="${IMMEDIATE_RETRY_DELAY_SECONDS:-60}"
NEXTFLOW_MODULE="${NEXTFLOW_MODULE:-}"
SAMTOOLS_MODULE="${SAMTOOLS_MODULE:-}"
STREAM_SMOKE_TEST="${STREAM_SMOKE_TEST:-0}"
PIPELINE_REPO_DIR="${PIPELINE_REPO_DIR:-}"

ACTIVE_STAGE=""
ACTIVE_CHILD_PID=""
ACTIVE_STAGE_LOG=""
REQUEUE_IN_PROGRESS=0
FAILED_STAGE="setup"
FAILURE_CLASS="UNKNOWN"
FAILURE_REASON=""
FAILURE_NONRETRYABLE=0
FIRST_FAILURE_EPOCH=""

log() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" >&2; }
die() { echo "ERROR: $*" >&2; exit 1; }
require_value() { [[ -n "${!1:-}" ]] || die "$1 must be explicitly supplied by the manager"; }
atomic_write() { local target=$1 tmp; mkdir -p "$(dirname "$target")"; tmp=$(mktemp "${target}.tmp.XXXXXX"); cat >"$tmp"; mv -f -- "$tmp" "$target"; }
truthy() { [[ "$1" == 1 || "$1" == true ]]; }

load_required_tools() {
    local tool module_name executable version
    for tool in nextflow samtools; do
        if [[ "$tool" == nextflow ]]; then module_name=$NEXTFLOW_MODULE; else module_name=$SAMTOOLS_MODULE; fi
        if [[ -n "$module_name" ]]; then
            command -v module >/dev/null 2>&1 || die "${tool^^}_MODULE is configured, but the module command is unavailable"
            module load "$module_name" || die "Unable to load module $module_name for $tool"
            hash -r
        elif ! command -v "$tool" >/dev/null 2>&1; then
            die "$tool is not in PATH and ${tool^^}_MODULE is empty"
        fi
        command -v "$tool" >/dev/null 2>&1 || die "$tool is unavailable after environment setup"
        executable=$(command -v "$tool")
        if [[ "$tool" == nextflow ]]; then
            version=$(nextflow -version 2>&1 | awk 'NF && !found {print; found=1}')
        else
            version=$(samtools --version 2>&1 | awk 'NF && !found {print; found=1}')
        fi
        [[ -n "$version" ]] || die "Unable to determine $tool version"
        log "INFO: ${tool}_executable=${executable}"
        log "INFO: ${tool}_version=${version}"
    done
}

if [[ -n "$PIPELINE_REPO_DIR" ]]; then
    REPO_DIR=$(realpath -m "$PIPELINE_REPO_DIR")
else
    REPO_DIR=$(realpath -m "$(dirname "${BASH_SOURCE[0]}")")
fi
case "$REPO_DIR" in
    /var/spool/slurmd/job*) die "Refusing to use Slurm spool directory as pipeline repository: $REPO_DIR" ;;
esac
export PIPELINE_REPO_DIR="$REPO_DIR"
log "INFO: pipeline_repo_dir=${REPO_DIR}"

required_pipeline_scripts=(
    "${REPO_DIR}/preprocessing.nf"
    "${REPO_DIR}/numt_detection/numt_end2end.nf"
    "${REPO_DIR}/primate_pipeline_numt_decoy_round1.nf"
    "${REPO_DIR}/primate_pipeline_round2_consensus_NUMT.nf"
)
for pipeline_script in "${required_pipeline_scripts[@]}"; do
    [[ -s "$pipeline_script" ]] || die "Required pipeline script is missing: $pipeline_script"
done

for variable in FULL_SAMPLE_LIST PRE_OUTPUT_DIR ROUND_OUTPUT_DIR NF_BASE_WORK_DIR NF_CONFIG_FILE GLOBAL_REF_DIR REF_DIR NUCLEAR_ONLY_REF_DIR; do require_value "$variable"; done
[[ -f "$FULL_SAMPLE_LIST" ]] || die "FULL_SAMPLE_LIST does not exist: $FULL_SAMPLE_LIST"
for directory in GLOBAL_REF_DIR REF_DIR NUCLEAR_ONLY_REF_DIR; do [[ -d "${!directory}" ]] || die "$directory is not a directory: ${!directory}"; done
NF_CONFIG_PATH=$([[ "$NF_CONFIG_FILE" = /* ]] && printf %s "$NF_CONFIG_FILE" || printf %s "${REPO_DIR}/${NF_CONFIG_FILE}")
[[ -s "$NF_CONFIG_PATH" ]] || die "NF_CONFIG_FILE does not exist or is empty: $NF_CONFIG_PATH"
mkdir -p "$NF_BASE_WORK_DIR" "${NF_BASE_WORK_DIR}/.sample_state" "${NF_BASE_WORK_DIR}/.locks" "${NF_BASE_WORK_DIR}/.manifests"
load_required_tools

normalize_manifest() {
    local source=$1 target=$2 tmp
    tmp="${target}.tmp.$$"
    awk -F '\t' 'BEGIN{OFS="\t"}
      {sub(/\r$/,"")} /^[[:space:]]*($|#)/{next}
      tolower($1)=="sample" || tolower($1)=="sample_id" {next}
      NF < 2 || $1=="" || $2=="" {next}
      { if ($1 in ref && ref[$1] != $2) {printf "conflicting reference for sample %s: %s versus %s\n",$1,ref[$1],$2 > "/dev/stderr"; bad=1; next}
        if (!($1 in ref)) {ref[$1]=$2; order[++n]=$1} }
      END {if (bad) exit 42; print "sample_id","reference_name"; for(i=1;i<=n;i++) print order[i],ref[order[i]]}' "$source" >"$tmp" || { rm -f "$tmp"; return 1; }
    [[ $(wc -l <"$tmp") -gt 1 ]] || { rm -f "$tmp"; die "No valid sample rows in $source"; }
    mv -f -- "$tmp" "$target"
}

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    submission_id="${SLURM_JOB_ID:-$(date -u +%Y%m%dT%H%M%SZ).$$}"
    NORMALIZED_SAMPLE_LIST="${NF_BASE_WORK_DIR}/.manifests/${submission_id}.samples.tsv"
    normalize_manifest "$FULL_SAMPLE_LIST" "$NORMALIZED_SAMPLE_LIST" || die "Unable to normalize sample manifest"
    NUM_SAMPLES=$(( $(wc -l <"$NORMALIZED_SAMPLE_LIST") - 1 ))
    if truthy "$STREAM_SMOKE_TEST"; then
        (( NUM_SAMPLES == 1 )) || die "STREAM_SMOKE_TEST=1 requires exactly one normalized sample (found ${NUM_SAMPLES})"
        MAX_CONCURRENT=1
        log "WARNING: STREAM_SMOKE_TEST=1; restricting this submission to one sample and one concurrent worker"
    fi
    export NORMALIZED_SAMPLE_LIST
    sbatch_args=(--export=ALL "--array=1-${NUM_SAMPLES}%${MAX_CONCURRENT}")
    [[ -z "${STREAM_PARTITION}" ]] || sbatch_args+=(--partition="${STREAM_PARTITION}")
    log "Submitting ${NUM_SAMPLES} sample workers from immutable manifest ${NORMALIZED_SAMPLE_LIST}"
    if ! submission=$(sbatch "${sbatch_args[@]}" "$0"); then
        die "sbatch failed while submitting the sample array"
    fi
    if [[ "$submission" =~ ^Submitted[[:space:]]+batch[[:space:]]+job[[:space:]]+([0-9]+)[[:space:]]*$ ]]; then
        submitted_job_id=${BASH_REMATCH[1]}
    else
        die "Unable to parse Slurm job ID from sbatch output: $submission"
    fi
    printf '%s\n' "$submission"
    log "INFO: submitted_array_job_id=${submitted_job_id}"
    log "INFO: normalized_manifest=${NORMALIZED_SAMPLE_LIST}"
    log "INFO: sample_count=${NUM_SAMPLES}"
    log "INFO: max_concurrent=${MAX_CONCURRENT}"
    exit 0
fi

[[ -n "${NORMALIZED_SAMPLE_LIST:-}" && -s "$NORMALIZED_SAMPLE_LIST" ]] || die "NORMALIZED_SAMPLE_LIST was not exported to the array worker"
SAMPLE_LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$NORMALIZED_SAMPLE_LIST")
IFS=$'\t' read -r SAMPLE_ID REF_NAME _ <<<"$SAMPLE_LINE"
if [[ -z "${SAMPLE_ID:-}" || ! "$SAMPLE_ID" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ || "$SAMPLE_ID" == *..* ]]; then
    die "UNSAFE_PATH: invalid SAMPLE_ID: ${SAMPLE_ID:-<empty>}"
fi
[[ -n "${REF_NAME:-}" ]] || die "MALFORMED_METADATA: missing reference_name"

SAMPLE_WORK_ROOT="${NF_BASE_WORK_DIR}/${SAMPLE_ID}"
root_real=$(realpath -m "$NF_BASE_WORK_DIR")
sample_real=$(realpath -m "$SAMPLE_WORK_ROOT")
case "$sample_real" in "$root_real"/*) ;; *) die "Unsafe sample work path";; esac
[[ "$sample_real" != "$root_real" && "$sample_real" != / ]] || die "Unsafe sample work path"

# A lock file is merely an inode; only flock determines whether a worker is live.
exec 9>"${NF_BASE_WORK_DIR}/.locks/${SAMPLE_ID}.lock"
flock -n 9 || die "Another worker is already using ${SAMPLE_ID}"
LOCK_HELD_SAMPLE="$SAMPLE_ID"

STATE_DIR="${NF_BASE_WORK_DIR}/.sample_state"
SAMPLE_TSV="${SAMPLE_WORK_ROOT}/inputs/${SAMPLE_ID}.sample.tsv"
METADATA_DIR="${SAMPLE_WORK_ROOT}/metadata"
LOG_DIR="${SAMPLE_WORK_ROOT}/logs"
PRE_WORK_DIR="${SAMPLE_WORK_ROOT}/pre"; NUMT_WORK_DIR="${SAMPLE_WORK_ROOT}/numt"
ROUND1_WORK_DIR="${SAMPLE_WORK_ROOT}/round1"; ROUND2_WORK_DIR="${SAMPLE_WORK_ROOT}/round2"
PRE_LAUNCH_DIR="${SAMPLE_WORK_ROOT}/nextflow_launch/pre"; NUMT_LAUNCH_DIR="${SAMPLE_WORK_ROOT}/nextflow_launch/numt"
ROUND1_LAUNCH_DIR="${SAMPLE_WORK_ROOT}/nextflow_launch/round1"; ROUND2_LAUNCH_DIR="${SAMPLE_WORK_ROOT}/nextflow_launch/round2"

safe_remove_sample_stage_work() {
    local id=$1 requested=$2 expected target
    [[ "$id" == "$SAMPLE_ID" && "${LOCK_HELD_SAMPLE:-}" == "$id" ]] || { log "refusing cleanup without the sample lock"; return 1; }
    target=$(realpath -m "$requested")
    case "$(basename "$target")" in pre|numt|round1|round2) ;; *) log "refusing unsafe sample stage cleanup: $target"; return 1;; esac
    expected=$(realpath -m "${NF_BASE_WORK_DIR}/${id}/$(basename "$target")")
    [[ "$target" == "$expected" && "$target" == "$root_real"/* && "$target" != / ]] || { log "refusing unsafe sample stage cleanup: $target"; return 1; }
    rm -rf -- "$target"
}
safe_remove_sample_work() {
    local id=$1 target expected
    [[ "$id" == "$SAMPLE_ID" && "${LOCK_HELD_SAMPLE:-}" == "$id" ]] || return 1
    target=$(realpath -m "${NF_BASE_WORK_DIR}/${id}"); expected=$(realpath -m "$SAMPLE_WORK_ROOT")
    [[ "$target" == "$expected" && "$target" == "$root_real"/* && "$target" != / && "$target" != "$root_real" ]] || { log "refusing unsafe sample cleanup: $target"; return 1; }
    case "$(basename "$target")" in .sample_state|.locks|.manifests) return 1;; esac
    rm -rf -- "$target"
}
clean_stage() { truthy "$CLEAN_VALIDATED_STAGE_WORK" && safe_remove_sample_stage_work "$SAMPLE_ID" "$1" || true; }

WGS_REF="${GLOBAL_REF_DIR}/${REF_NAME}.fa"; CHRM_REF="${REF_DIR}/${REF_NAME}.fa"
NUCLEAR_REF="${NUCLEAR_ONLY_REF_DIR}/${REF_NAME}.nuclear_only.fa"
[[ -s "$NUCLEAR_REF" ]] || NUCLEAR_REF="${NUCLEAR_ONLY_REF_DIR}/${REF_NAME}.fa"
for ref in "$WGS_REF" "$CHRM_REF" "$NUCLEAR_REF"; do
    if [[ ! -s "$ref" ]]; then FAILURE_CLASS=MISSING_REFERENCE; FAILURE_REASON="Missing reference: $ref"; FAILURE_NONRETRYABLE=1; break; fi
done

signature() { local p=$1; if [[ -e "$p" ]]; then printf '%s\t%s\t%s\n' "$(realpath -m "$p")" "$(stat -c %s "$p")" "$(stat -c %Y "$p")"; else printf '%s\tMISSING\n' "$(realpath -m "$p")"; fi; }
build_fingerprint() {
    local out=$1 manifest_sha config_sha commit reference_fp parameters_fp combined
    manifest_sha=$(printf '%s\n' "$SAMPLE_LINE" | sha256sum | awk '{print $1}') || return 1
    config_sha=$(sha256sum "$NF_CONFIG_PATH" | awk '{print $1}') || return 1
    commit=$(git -C "$REPO_DIR" rev-parse HEAD 2>/dev/null || printf unknown)
    reference_fp=$({ for r in "$WGS_REF" "$CHRM_REF" "$NUCLEAR_REF"; do signature "$r"; for x in "$r.fai" "$r.dict" "${r%.*}.dict" "$r.amb" "$r.ann" "$r.bwt" "$r.pac" "$r.sa"; do signature "$x"; done; done; } | sha256sum | awk '{print $1}') || return 1
    parameters_fp=$(printf '%s\n' "ENABLE_CHUNKED_ALIGNMENT=${ENABLE_CHUNKED_ALIGNMENT:-}" "FASTQ_SIZE_THRESHOLD_GB=${FASTQ_SIZE_THRESHOLD_GB:-}" "READS_PER_CHUNK=${READS_PER_CHUNK:-}" "NF_CONFIG_FILE=$NF_CONFIG_PATH" "pipeline_profile=$PIPELINE_PROFILE" "GLOBAL_REF_DIR=$(realpath -m "$GLOBAL_REF_DIR")" "REF_DIR=$(realpath -m "$REF_DIR")" "NUCLEAR_ONLY_REF_DIR=$(realpath -m "$NUCLEAR_ONLY_REF_DIR")" | sha256sum | awk '{print $1}') || return 1
    combined=$(printf '%s\n' "$manifest_sha" "$config_sha" "$commit" "$reference_fp" "$parameters_fp" | sha256sum | awk '{print $1}') || return 1
    { printf 'key\tvalue\n'; printf 'sample_manifest_sha256\t%s\npipeline_config_sha256\t%s\npipeline_git_commit\t%s\nreference_fingerprint\t%s\nimportant_parameters_sha256\t%s\ncombined_fingerprint\t%s\n' "$manifest_sha" "$config_sha" "$commit" "$reference_fp" "$parameters_fp" "$combined"; } >"$out"
}

fingerprint_tmp=$(mktemp "${NF_BASE_WORK_DIR}/.fingerprint.${SAMPLE_ID}.XXXXXX")
if ! build_fingerprint "$fingerprint_tmp"; then FAILURE_CLASS=FINGERPRINT_GENERATION_FAILED; FAILURE_REASON="Fingerprint generation failed"; FAILURE_NONRETRYABLE=1; fi
fingerprint_match=0; resume_mode=fresh
if [[ -s "$fingerprint_tmp" && -f "${METADATA_DIR}/fingerprint.tsv" ]] && cmp -s "$fingerprint_tmp" "${METADATA_DIR}/fingerprint.tsv"; then fingerprint_match=1; resume_mode=resume
elif [[ -e "$SAMPLE_WORK_ROOT" ]]; then
    stale="${NF_BASE_WORK_DIR}/${SAMPLE_ID}.stale.$(date -u +%Y%m%dT%H%M%SZ).$$"
    mv -- "$SAMPLE_WORK_ROOT" "$stale"
    log "Fingerprint changed; archived old cache at $stale"
fi
mkdir -p "$METADATA_DIR" "$LOG_DIR" "$(dirname "$SAMPLE_TSV")" "$PRE_LAUNCH_DIR" "$NUMT_LAUNCH_DIR" "$ROUND1_LAUNCH_DIR" "$ROUND2_LAUNCH_DIR"
[[ -s "$fingerprint_tmp" ]] && mv -f "$fingerprint_tmp" "${METADATA_DIR}/fingerprint.tsv" || rm -f "$fingerprint_tmp"
printf '%s\n' "$SAMPLE_LINE" >"$SAMPLE_TSV"
COMBINED_FINGERPRINT=$(awk -F '\t' '$1=="combined_fingerprint"{print $2}' "${METADATA_DIR}/fingerprint.tsv" 2>/dev/null || true)
# Preserve queue age across later manager submissions; cleanup policy sorts by
# this value rather than mutable directory mtimes or the most recent attempt.
if [[ -s "${STATE_DIR}/${SAMPLE_ID}.failure.tsv" ]]; then
    FIRST_FAILURE_EPOCH=$(awk -F '\t' '$1=="first_failure_epoch"{print $2; exit}' "${STATE_DIR}/${SAMPLE_ID}.failure.tsv")
fi
log "INFO: sample_fingerprint_match=${fingerprint_match}"
log "INFO: resume_mode=${resume_mode}"
log "INFO: sample_work_root=${SAMPLE_WORK_ROOT}"

current_slurm_element_id() {
    if [[ -n "${SLURM_ARRAY_JOB_ID:-}" && -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        printf '%s_%s\n' "${SLURM_ARRAY_JOB_ID}" "${SLURM_ARRAY_TASK_ID}"
    else
        printf '%s\n' "${SLURM_JOB_ID:-}"
    fi
}
write_running_marker() {
    local worker_state=$1 current_attempt=$2
    atomic_write "${STATE_DIR}/${SAMPLE_ID}.running.tsv" <<EOF
sample_id\t${SAMPLE_ID}
reference_name\t${REF_NAME}
work_root\t${SAMPLE_WORK_ROOT}
worker_state\t${worker_state}
active_stage\t${ACTIVE_STAGE:-none}
attempt\t${current_attempt}
slurm_job_id\t${SLURM_JOB_ID:-}
array_job_id\t${SLURM_ARRAY_JOB_ID:-}
array_task_id\t${SLURM_ARRAY_TASK_ID:-}
updated_epoch\t$(date +%s)
updated_time\t$(date -u +%FT%TZ)
EOF
}
handle_sample_requeue() {
    [[ "$REQUEUE_IN_PROGRESS" == 0 ]] || return 0; REQUEUE_IN_PROGRESS=1
    log "INFO: Received USR1 before walltime"
    log "INFO: sample_id=${SAMPLE_ID}"
    log "INFO: active_stage=${ACTIVE_STAGE:-none}"
    log "INFO: preserving_sample_work_root=${SAMPLE_WORK_ROOT}"
    write_running_marker REQUEUE_REQUESTED "${attempt:-0}"
    atomic_write "${STATE_DIR}/${SAMPLE_ID}.requeue.tsv" <<EOF
sample_id\t${SAMPLE_ID}
active_stage\t${ACTIVE_STAGE:-none}
reason\tTIMEOUT_SIGNAL
resume_eligible\t1
work_root\t${SAMPLE_WORK_ROOT}
slurm_job_id\t${SLURM_JOB_ID:-}
array_job_id\t${SLURM_ARRAY_JOB_ID:-}
array_task_id\t${SLURM_ARRAY_TASK_ID:-}
time\t$(date -u +%FT%TZ)
EOF
    if [[ -n "$ACTIVE_CHILD_PID" ]] && kill -0 "$ACTIVE_CHILD_PID" 2>/dev/null; then
        kill -TERM "$ACTIVE_CHILD_PID" 2>/dev/null || true
        for _ in $(seq 1 60); do kill -0 "$ACTIVE_CHILD_PID" 2>/dev/null || break; sleep 2; done
        if kill -0 "$ACTIVE_CHILD_PID" 2>/dev/null; then kill -KILL "$ACTIVE_CHILD_PID" 2>/dev/null || true; fi
        wait "$ACTIVE_CHILD_PID" 2>/dev/null || true
        ACTIVE_CHILD_PID=""
    fi
    element_id=$(current_slurm_element_id)
    if [[ -z "$element_id" ]]; then
        log "ERROR: cannot requeue: current Slurm element ID is empty"
        exit 75
    fi
    log "INFO: requeueing_array_element=${element_id}"
    if ! scontrol requeue "$element_id"; then
        log "ERROR: scontrol requeue failed for ${element_id}"
        exit 75
    fi
    exit 0
}
trap handle_sample_requeue USR1

run_stage() {
    local stage=$1 stage_log=$2; shift 2
    ACTIVE_STAGE=$stage; ACTIVE_STAGE_LOG=$stage_log
    mkdir -p "$(dirname "$stage_log")"
    set +e
    ( cd "${SAMPLE_WORK_ROOT}/nextflow_launch/${stage}"; "$@" ) > >(tee -a "$stage_log") 2> >(tee -a "$stage_log" >&2) &
    ACTIVE_CHILD_PID=$!; wait "$ACTIVE_CHILD_PID"; local rc=$?; ACTIVE_CHILD_PID=""; set -e
    (( REQUEUE_IN_PROGRESS == 0 )) || return 125
    if (( rc != 0 )); then
        FAILED_STAGE=$stage
        # DOWNLOAD_FASTQ emits this stable marker for deterministic ENA
        # metadata/layout failures.  Preserve it instead of replacing it with
        # a generic CRAM validation error, and suppress the outer retry.
        if [[ "$stage" == pre && -f "$stage_log" ]]; then
            local root_cause
            root_cause=$(awk '/DETERMINISTIC_FASTQ_FAILURE class=/{line=$0} END{print line}' "$stage_log")
            if [[ "$root_cause" =~ class=([^[:space:]]+)[[:space:]]+run=([^[:space:]]+)[[:space:]]+reason=(.*)$ ]]; then
                FAILURE_CLASS=${BASH_REMATCH[1]}
                FAILURE_REASON="${BASH_REMATCH[2]} ${BASH_REMATCH[3]}"
                FAILURE_NONRETRYABLE=1
            fi
        fi
        if [[ -z "$FAILURE_CLASS" || "$FAILURE_CLASS" == UNKNOWN ]]; then
            FAILURE_CLASS=STAGE_FAILED; FAILURE_REASON="${stage} Nextflow exited ${rc}"
        fi
        return "$rc"
    fi
}
nf() { nextflow -C "$NF_CONFIG_PATH" run "$1" -profile "$PIPELINE_PROFILE" -resume -w "$2" "${@:3}"; }
file_nonempty() { [[ -s "$1" ]]; }
find_nonempty() { find "$1" -type f -name "$2" -size +0c -print -quit 2>/dev/null | grep -q .; }

CRAM_PATH="${PRE_OUTPUT_DIR}/${SAMPLE_ID}/alignment/${SAMPLE_ID}.cram"; CRAI_PATH="${CRAM_PATH}.crai"
NUMT_BED="${NUMT_BESTHIT_OUTDIR}/${SAMPLE_ID}.highconf_numt.bed"
ROUND1_BAM="${ROUND1_OUTDIR}/${SAMPLE_ID}/round_1/candidate_reads/${SAMPLE_ID}.with_mates.bam"
ROUND1_VCF_DIR="${ROUND1_OUTDIR}/${SAMPLE_ID}/round_1_variant_calling_decoy"
ROUND1_NUMT_FA="${ROUND1_OUTDIR}/${SAMPLE_ID}/round_1/numt_decoy_ref/${SAMPLE_ID}.original_numt.fa"
ROUND1_NUMT_VCF="${ROUND1_OUTDIR}/${SAMPLE_ID}/round_1/numt_decoy_variant_calling/${SAMPLE_ID}.numt_decoy.raw.vcf.gz"
ROUND2_VCF_DIR="${ROUND_OUTPUT_DIR}/${SAMPLE_ID}/round_2_variant_calling_original_coords"
ROUND2_VCF="${ROUND2_VCF_DIR}/${SAMPLE_ID}.round2.original_coords.clean.final.split.vcf.gz"
ROUND2_COVERAGE="${ROUND2_VCF_DIR}/${SAMPLE_ID}.round2.original_coords.per_base_coverage.tsv"
ROUND2_MTCN="${ROUND_OUTPUT_DIR}/${SAMPLE_ID}/round_2/mtcn/${SAMPLE_ID}.round2.mtcn.tsv"

pre_complete() { [[ -s "$CRAM_PATH" && -s "$CRAI_PATH" ]] && samtools quickcheck "$CRAM_PATH"; }
numt_complete() { [[ -e "$NUMT_BED" ]]; }
round1_complete() { [[ -s "$ROUND1_BAM" && -e "$ROUND1_NUMT_FA" && -s "$ROUND1_NUMT_VCF" && -s "${ROUND1_NUMT_VCF}.tbi" ]] && find_nonempty "$ROUND1_VCF_DIR" "${SAMPLE_ID}.numt_decoy.clean.final.split.vcf"; }
round2_complete() { [[ -s "$ROUND2_VCF" && -s "$ROUND2_COVERAGE" && -s "$ROUND2_MTCN" ]]; }
numt_decoy_coverage_path() { find "$ROUND1_VCF_DIR" -type f \( -name '*per_base_coverage*' -o -name '*per-base*coverage*' \) -size +0c -print -quit 2>/dev/null; }
final_outputs_complete() {
    [[ -s "$CRAM_PATH" && -s "$CRAI_PATH" ]] || return 1
    samtools quickcheck "$CRAM_PATH" || return 1
    [[ -s "$ROUND2_VCF" ]] && gzip -t "$ROUND2_VCF" || return 1
    [[ -s "$ROUND2_COVERAGE" ]] || return 1
    [[ -n "$(numt_decoy_coverage_path)" ]] || return 1
    [[ -s "$ROUND2_MTCN" ]] || return 1
}

write_numt_config() {
    [[ -s "${CHRM_REF}.fai" ]] || { FAILURE_CLASS=MISSING_REFERENCE; FAILURE_REASON="Missing ${CHRM_REF}.fai"; FAILURE_NONRETRYABLE=1; return 1; }
    local mt_contig mt_length config="${METADATA_DIR}/${SAMPLE_ID}.numt.config"
    read -r mt_contig mt_length _ <"${CHRM_REF}.fai"
    cat >"$config" <<EOF
SAMPLE=${SAMPLE_ID}
SAMPLE_ID=${SAMPLE_ID}
SPECIES_NAME=${REF_NAME}
REF_NAME=${REF_NAME}
SAMPLES_TSV=${SAMPLE_TSV}
INPUT_BAM_CRAM=${CRAM_PATH}
INPUT_INDEX=${CRAI_PATH}
INPUT_BAM_CRAM_ALT=${CRAM_PATH}
INPUT_INDEX_ALT=${CRAI_PATH}
MT_CONTIG=${mt_contig}
MT_LENGTH=${mt_length}
WGS_REF=${WGS_REF}
NUCLEAR_REF=${NUCLEAR_REF}
CHRM_REF=${CHRM_REF}
CRAM_ROOT_1=${PRE_OUTPUT_DIR}
CRAM_ROOT_2=${PRE_OUTPUT_DIR}
WHOLE_REF_DIR=${GLOBAL_REF_DIR}
NUCLEAR_ONLY_REF_DIR=${NUCLEAR_ONLY_REF_DIR}
CHRM_REF_DIR=${REF_DIR}
OUTDIR=${NUMT_DISCOVERY_OUTROOT}
DISCOVERY_OUTDIR=${NUMT_DISCOVERY_OUTROOT}
DISCOVERY_OUTROOT=${NUMT_DISCOVERY_OUTROOT}
BESTHIT_OUTDIR=${NUMT_BESTHIT_OUTDIR}
EOF
    NUMT_CONFIG=$config
}

run_sample_chain() {
    if (( FAILURE_NONRETRYABLE )); then FAILED_STAGE=setup; return 1; fi
    mkdir -p "$PRE_WORK_DIR" "$NUMT_WORK_DIR" "$ROUND1_WORK_DIR" "$ROUND2_WORK_DIR" "$NUMT_DISCOVERY_OUTROOT" "$NUMT_BESTHIT_OUTDIR"
    if ! pre_complete; then run_stage pre "${LOG_DIR}/pre.log" nf "${REPO_DIR}/preprocessing.nf" "$PRE_WORK_DIR" --sample_tsv "$SAMPLE_TSV" --outdir "$PRE_OUTPUT_DIR" || return; fi
    pre_complete || { FAILED_STAGE=pre; FAILURE_CLASS=OUTPUT_INCOMPLETE; FAILURE_REASON="CRAM/CRAI validation failed"; return 1; }; clean_stage "$PRE_WORK_DIR"
    if ! numt_complete; then write_numt_config || return 1; run_stage numt "${LOG_DIR}/numt.log" nf "${REPO_DIR}/numt_detection/numt_end2end.nf" "$NUMT_WORK_DIR" --numt_config "$NUMT_CONFIG" || return; fi
    numt_complete || { FAILED_STAGE=numt; FAILURE_CLASS=OUTPUT_INCOMPLETE; FAILURE_REASON="NUMT BED missing"; return 1; }; clean_stage "$NUMT_WORK_DIR"
    if ! round1_complete; then run_stage round1 "${LOG_DIR}/round1.log" nf "${REPO_DIR}/primate_pipeline_numt_decoy_round1.nf" "$ROUND1_WORK_DIR" --sample_tsv "$SAMPLE_TSV" --outdir "$ROUND1_OUTDIR" --cram_dirs "$PRE_OUTPUT_DIR" --numt_bed_dir "$NUMT_BESTHIT_OUTDIR" || return; fi
    round1_complete || { FAILED_STAGE=round1; FAILURE_CLASS=OUTPUT_INCOMPLETE; FAILURE_REASON="Round 1 outputs incomplete"; return 1; }; clean_stage "$ROUND1_WORK_DIR"
    if ! round2_complete; then run_stage round2 "${LOG_DIR}/round2.log" nf "${REPO_DIR}/primate_pipeline_round2_consensus_NUMT.nf" "$ROUND2_WORK_DIR" --sample_tsv "$SAMPLE_TSV" --outdir "$ROUND_OUTPUT_DIR" --round1_outdir "$ROUND1_OUTDIR" || return; fi
    # Do not clean round2 until the complete eight-part final validation succeeds.
    final_outputs_complete || { FAILED_STAGE=round2; FAILURE_CLASS=OUTPUT_INCOMPLETE; FAILURE_REASON="Final eight-part output validation failed"; return 1; }
}

record_attempt() { printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$(date +%s)" "$FAILED_STAGE" "$FAILURE_CLASS" "$FAILURE_REASON" >>"${METADATA_DIR}/attempts.tsv"; }
classify_sample_failure() { case "$FAILURE_CLASS" in MISSING_REFERENCE|MALFORMED_METADATA|UNSUPPORTED_REFERENCE|UNSAFE_PATH|FINGERPRINT_GENERATION_FAILED|UNSUPPORTED_FASTQ_LAYOUT|AMBIGUOUS_FASTQ_LAYOUT|MD5_URL_COUNT_MISMATCH|MALFORMED_ENA_METADATA|MISSING_R1|MISSING_R2|DUPLICATE_R1|DUPLICATE_R2) FAILURE_NONRETRYABLE=1;; esac; }
collect_diagnostics() {
    local d="${METADATA_DIR}/diagnostics.$(date -u +%Y%m%dT%H%M%SZ)"; mkdir -p "$d"
    [[ -f "$ACTIVE_STAGE_LOG" ]] && tail -n 500 "$ACTIVE_STAGE_LOG" >"$d/stage.log.tail" || true
    [[ -f "${SAMPLE_WORK_ROOT}/nextflow_launch/${FAILED_STAGE}/.nextflow.log" ]] && tail -n 500 "${SAMPLE_WORK_ROOT}/nextflow_launch/${FAILED_STAGE}/.nextflow.log" >"$d/nextflow.log.tail" || true
    { scontrol show job "${SLURM_JOB_ID:-}" 2>&1 || true; } >"$d/slurm.txt"
    df -h >"$d/df-h.txt" 2>&1 || true; df -ih >"$d/df-ih.txt" 2>&1 || true; du -sb "$SAMPLE_WORK_ROOT" >"$d/du-sb.txt" 2>&1 || true
    cp "${METADATA_DIR}/fingerprint.tsv" "$d/" 2>/dev/null || true; printf %s "$d"
}
write_failure_marker() {
    local now=$1 diagnostics=$2 work_bytes
    work_bytes=$(du -sb "$SAMPLE_WORK_ROOT" 2>/dev/null | awk '{print $1}'); work_bytes=${work_bytes:-0}
    atomic_write "${STATE_DIR}/${SAMPLE_ID}.failure.tsv" <<EOF
sample_id\t${SAMPLE_ID}
reference_name\t${REF_NAME}
worker_state\tTERMINAL_FAILED
failed_stage\t${FAILED_STAGE}
failure_class\t${FAILURE_CLASS}
failure_reason\t${FAILURE_REASON//$'\n'/ }
attempt_count\t${attempt}
first_failure_time\t$(date -u -d "@${FIRST_FAILURE_EPOCH}" +%FT%TZ)
last_failure_time\t$(date -u -d "@${now}" +%FT%TZ)
first_failure_epoch\t${FIRST_FAILURE_EPOCH}
last_failure_epoch\t${now}
last_attempt_epoch\t${now}
work_root\t${SAMPLE_WORK_ROOT}
work_bytes\t${work_bytes}
fingerprint\t${COMBINED_FINGERPRINT}
resume_possible\t1
diagnostics_path\t${diagnostics}
slurm_job_id\t${SLURM_JOB_ID:-}
array_job_id\t${SLURM_ARRAY_JOB_ID:-}
array_task_id\t${SLURM_ARRAY_TASK_ID:-}
EOF
}

attempt=0; max_attempts=$((1 + IMMEDIATE_SAMPLE_RETRIES)); success=0
write_running_marker RUNNING "$attempt"
while (( attempt < max_attempts )); do
    attempt=$((attempt + 1))
    if run_sample_chain; then success=1; break; fi
    (( REQUEUE_IN_PROGRESS == 0 )) || exit 0
    now=$(date +%s); [[ -n "$FIRST_FAILURE_EPOCH" ]] || FIRST_FAILURE_EPOCH=$now
    record_attempt "$attempt"; classify_sample_failure
    (( FAILURE_NONRETRYABLE == 0 )) || break
    if (( attempt < max_attempts )); then
        write_running_marker IMMEDIATE_RETRY "$attempt"
        log "Immediate retry ${attempt}/${IMMEDIATE_SAMPLE_RETRIES}"
        sleep "$IMMEDIATE_RETRY_DELAY_SECONDS"
    fi
done

if (( success )); then
    # Persistence ordering is deliberate: state survives optional work-root removal.
    atomic_write "${STATE_DIR}/${SAMPLE_ID}.complete.tsv" <<EOF
sample_id\t${SAMPLE_ID}
reference_name\t${REF_NAME}
fingerprint\t${COMBINED_FINGERPRINT}
completed_epoch\t$(date +%s)
completed_time\t$(date -u +%FT%TZ)
EOF
    rm -f -- "${STATE_DIR}/${SAMPLE_ID}.failure.tsv" "${STATE_DIR}/${SAMPLE_ID}.requeue.tsv" "${STATE_DIR}/${SAMPLE_ID}.running.tsv"
    clean_stage "$PRE_WORK_DIR"; clean_stage "$NUMT_WORK_DIR"; clean_stage "$ROUND1_WORK_DIR"; clean_stage "$ROUND2_WORK_DIR"
    truthy "$REMOVE_SAMPLE_ROOT_ON_SUCCESS" && safe_remove_sample_work "$SAMPLE_ID"
    log "Completed per-sample streaming pipeline for ${SAMPLE_ID}"
    exit 0
fi

now=$(date +%s); [[ -n "$FIRST_FAILURE_EPOCH" ]] || FIRST_FAILURE_EPOCH=$now
diagnostics=$(collect_diagnostics)
write_failure_marker "$now" "$diagnostics"
rm -f -- "${STATE_DIR}/${SAMPLE_ID}.running.tsv"
log "Sample ${SAMPLE_ID} failed after ${attempt} attempt(s); retaining ${FAILED_STAGE} work"
exit 1
