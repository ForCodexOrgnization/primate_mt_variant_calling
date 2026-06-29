#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: bash run_numt_one_from_3col.sh <samples.tsv> <task_id> <results_root> <nuclear_only_ref_dir> [threads] [min_mapq] [min_depth] [min_reads] [min_len] [merge_gap] [pad]" >&2
  exit 1
fi

SAMPLES_TSV="$1"
TASK_ID="$2"
RESULTS_ROOT="$3"
NUCLEAR_ONLY_REF_DIR="$4"

# Optional: directory for checking existing discovery results.
# If not provided via env, default to RESULTS_ROOT.
DISCOVERY_OUTROOT="${DISCOVERY_OUTROOT:-${RESULTS_ROOT}}"

THREADS="${5:-8}"
MIN_MAPQ="${6:-20}"
MIN_DEPTH="${7:-3}"
MIN_READS="${8:-5}"
MIN_LEN="${9:-100}"
MERGE_GAP="${10:-50}"
PAD="${11:-500}"

CRAM_ROOT_1="/nfs/roberts/pi/pi_njl27/lt692/primate_ref_files/nuclear_only_refs"
CRAM_ROOT_2="/home/lt692/scratch_pi_njl27/lt692/primate_results"

REF_ROOT="/nfs/roberts/pi/pi_njl27/lt692/primate_ref_files/whole_genomes_final"
MT_CONTIG="chrM"

LINE=$(sed -n "${TASK_ID}p" "${SAMPLES_TSV}")
if [[ -z "${LINE}" ]]; then
  echo "ERROR: No line found for task_id=${TASK_ID}" >&2
  exit 1
fi

PARSED_FIELDS=$(printf '%s\n' "$LINE" | awk -F'[\t,]' '
{
  n=0
  for (i=1; i<=NF; i++) {
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", $i)
    if ($i != "") {
      n++
      f[n]=$i
      if (n==3) break
    }
  }
  if (n < 3) {
    n=split($0, raw, /[[:space:]]+/)
    m=0
    for (j=1; j<=n; j++) {
      if (raw[j] != "") {
        m++
        f[m]=raw[j]
        if (m==3) break
      }
    }
    if (m >= 1) print f[1]
    if (m >= 2) print f[2]
    if (m >= 3) print f[3]
  } else {
    print f[1]
    print f[2]
    print f[3]
  }
}')
SAMPLE_ID=$(printf '%s\n' "$PARSED_FIELDS" | sed -n '1p')
SPECIES_NAME=$(printf '%s\n' "$PARSED_FIELDS" | sed -n '2p')
REF_NAME=$(printf '%s\n' "$PARSED_FIELDS" | sed -n '3p')

[[ -n "${SAMPLE_ID}" ]] || { echo "ERROR: sample_id empty at task ${TASK_ID}" >&2; exit 1; }
[[ -n "${REF_NAME}" ]] || { echo "ERROR: ref_name empty at task ${TASK_ID}" >&2; exit 1; }

resolve_cram_in_one_root() {
  local sample_id="$1"
  local cram_root="$2"

  local hit
  hit=$(find "${cram_root}" -type f \( -name "${sample_id}.cram" -o -name "${sample_id}*.cram" \) | head -n 1 || true)

  if [[ -n "${hit}" ]]; then
    echo "${hit}"
    return 0
  fi

  return 1
}

resolve_cram() {
  local sample_id="$1"

  local hit=""

  hit=$(resolve_cram_in_one_root "${sample_id}" "${CRAM_ROOT_1}" || true)
  if [[ -n "${hit}" ]]; then
    echo "${hit}"
    return 0
  fi

  hit=$(resolve_cram_in_one_root "${sample_id}" "${CRAM_ROOT_2}" || true)
  if [[ -n "${hit}" ]]; then
    echo "${hit}"
    return 0
  fi

  return 1
}

resolve_full_ref() {
  local ref_name="$1"
  local ref_root="$2"

  for ext in fa fasta fna; do
    if [[ -f "${ref_root}/${ref_name}.${ext}" ]]; then
      echo "${ref_root}/${ref_name}.${ext}"
      return 0
    fi
  done

  local hit
  hit=$(find "${ref_root}" -maxdepth 1 -type f \( -name "${ref_name}.fa" -o -name "${ref_name}.fasta" -o -name "${ref_name}.fna" -o -name "${ref_name}_genomic.fa" -o -name "${ref_name}_genomic.fna" -o -name "${ref_name}*.fa" -o -name "${ref_name}*.fasta" -o -name "${ref_name}*.fna" \) | head -n 1 || true)
  [[ -n "${hit}" ]] || return 1
  echo "${hit}"
}

resolve_nuclear_ref() {
  local ref_name="$1"
  local ref_root="$2"

  if [[ -f "${ref_root}/${ref_name}.nuclear_only.fa" ]]; then
    echo "${ref_root}/${ref_name}.nuclear_only.fa"
    return 0
  fi

  local hit
  hit=$(find "${ref_root}" -maxdepth 1 -type f -name "${ref_name}*.nuclear_only.fa" | head -n 1 || true)
  [[ -n "${hit}" ]] || return 1
  echo "${hit}"
}

INPUT_CRAM=$(resolve_cram "${SAMPLE_ID}") || {
  echo "ERROR: cannot find CRAM for ${SAMPLE_ID} under either:" >&2
  echo "  ${CRAM_ROOT_1}" >&2
  echo "  ${CRAM_ROOT_2}" >&2
  exit 1
}

if [[ -f "${INPUT_CRAM}.crai" ]]; then
  INPUT_CRAI="${INPUT_CRAM}.crai"
elif [[ -f "${INPUT_CRAM%.cram}.cram.crai" ]]; then
  INPUT_CRAI="${INPUT_CRAM%.cram}.cram.crai"
else
  echo "ERROR: cannot find CRAI for ${INPUT_CRAM}" >&2
  exit 1
fi

FULL_REF=$(resolve_full_ref "${REF_NAME}" "${REF_ROOT}") || {
  echo "ERROR: cannot find full reference fasta for ${REF_NAME}" >&2
  exit 1
}

NUCLEAR_REF=$(resolve_nuclear_ref "${REF_NAME}" "${NUCLEAR_ONLY_REF_DIR}") || {
  echo "ERROR: cannot find nuclear-only fasta for ${REF_NAME} under ${NUCLEAR_ONLY_REF_DIR}" >&2
  exit 1
}

samtools faidx "${FULL_REF}" >/dev/null 2>&1 || true

MT_LENGTH=$(awk -v mt="${MT_CONTIG}" '$1==mt {print $2}' "${FULL_REF}.fai" | head -n 1)
if [[ -z "${MT_LENGTH}" ]]; then
  echo "ERROR: cannot find ${MT_CONTIG} in ${FULL_REF}.fai" >&2
  exit 1
fi

OUTDIR="${RESULTS_ROOT}/${SAMPLE_ID}"
DISCOVERY_SAMPLE_OUTDIR="${DISCOVERY_OUTROOT}/${SAMPLE_ID}"

# Skip this sample when corresponding outputs already exist in DISCOVERY_OUTROOT.
if compgen -G "${DISCOVERY_SAMPLE_OUTDIR}"/*.numt_candidates.bed > /dev/null \
  || compgen -G "${DISCOVERY_SAMPLE_OUTDIR}"/*.numt_candidates.tsv > /dev/null; then
  echo "Skip ${SAMPLE_ID}: existing discovery outputs found in ${DISCOVERY_SAMPLE_OUTDIR}"
  exit 0
fi

echo "========================================"
echo "TASK_ID      : ${TASK_ID}"
echo "SAMPLE_ID    : ${SAMPLE_ID}"
echo "SPECIES      : ${SPECIES_NAME}"
echo "REF_NAME     : ${REF_NAME}"
echo "INPUT_CRAM   : ${INPUT_CRAM}"
echo "INPUT_CRAI   : ${INPUT_CRAI}"
echo "FULL_REF     : ${FULL_REF}"
echo "NUCLEAR_REF  : ${NUCLEAR_REF}"
echo "MT_CONTIG    : ${MT_CONTIG}"
echo "MT_LENGTH    : ${MT_LENGTH}"
echo "OUTDIR       : ${OUTDIR}"
echo "DISCOVERY_OUTROOT : ${DISCOVERY_OUTROOT}"
echo "========================================"

bash run_numt_discovery.sh \
  --input "${INPUT_CRAM}" \
  --index "${INPUT_CRAI}" \
  --sample "${SAMPLE_ID}" \
  --mt-contig "${MT_CONTIG}" \
  --mt-length "${MT_LENGTH}" \
  --wgs-ref "${FULL_REF}" \
  --nuclear-ref "${NUCLEAR_REF}" \
  --outdir "${OUTDIR}" \
  --threads "${THREADS}" \
  --min-mapq "${MIN_MAPQ}" \
  --min-depth "${MIN_DEPTH}" \
  --min-reads "${MIN_READS}" \
  --min-len "${MIN_LEN}" \
  --merge-gap "${MERGE_GAP}" \
  --pad "${PAD}"
