#!/usr/bin/env bash
#SBATCH --job-name=numt_e2e
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/numt_e2e_%A_%a.out
#SBATCH --error=logs/numt_e2e_%A_%a.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}"
source "${REPO_DIR}/load_numt_modules.sh"
cd "$REPO_DIR"

usage() {
  cat <<USAGE
Usage:
  bash submit_numt_end2end_array.sh --config numt_pipeline.config [--concurrent 50]

Required config keys for array mode:
  SAMPLES_TSV            # 3-column tsv: sample_id, real_species, ref_species
  CRAM_ROOT_1
  CRAM_ROOT_2
  WHOLE_REF_DIR
  NUCLEAR_ONLY_REF_DIR
  CHRM_REF_DIR
  DISCOVERY_OUTROOT      # parent dir for per-sample discovery output
  BESTHIT_OUTDIR
USAGE
}

CONFIG=""
CONCURRENT=50
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --concurrent) CONCURRENT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$CONFIG" ]] || { echo "ERROR: --config required" >&2; exit 1; }
[[ -s "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG"

: "${SAMPLES_TSV:=${SAMPLES_FILE:-}}"
: "${SAMPLES_TSV:?missing SAMPLES_TSV in config}"
: "${CRAM_ROOT_1:?missing CRAM_ROOT_1 in config}"
: "${CRAM_ROOT_2:?missing CRAM_ROOT_2 in config}"
: "${WHOLE_REF_DIR:?missing WHOLE_REF_DIR in config}"
: "${NUCLEAR_ONLY_REF_DIR:?missing NUCLEAR_ONLY_REF_DIR in config}"
: "${CHRM_REF_DIR:?missing CHRM_REF_DIR in config}"
: "${DISCOVERY_OUTROOT:?missing DISCOVERY_OUTROOT in config}"
: "${BESTHIT_OUTDIR:?missing BESTHIT_OUTDIR in config}"

[[ -s "$SAMPLES_TSV" ]] || { echo "ERROR: samples TSV not found: $SAMPLES_TSV" >&2; exit 1; }
mkdir -p logs "$DISCOVERY_OUTROOT" "$BESTHIT_OUTDIR"

N=$(wc -l < "$SAMPLES_TSV")
[[ "$N" -gt 0 ]] || { echo "ERROR: SAMPLES_TSV is empty: $SAMPLES_TSV" >&2; exit 1; }

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  sbatch --array=1-${N}%${CONCURRENT} "$0" --config "$CONFIG" --concurrent "$CONCURRENT"
  exit 0
fi

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_TSV" || true)
[[ -n "$LINE" ]] || { echo "ERROR: empty line for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
# Support TSV/CSV/whitespace-delimited sample sheets to keep single-sample and Slurm modes consistent.
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
REAL_SPECIES=$(printf '%s\n' "$PARSED_FIELDS" | sed -n '2p')
REF_SPECIES=$(printf '%s\n' "$PARSED_FIELDS" | sed -n '3p')
[[ -n "$SAMPLE_ID" ]] || { echo "ERROR: empty sample at task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
[[ -n "$REF_SPECIES" ]] || { echo "ERROR: empty ref species at task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }

resolve_cram_in_root() {
  local sid="$1"; local root="$2"
  find "$root" -type f \( -name "${sid}.cram" -o -name "${sid}*.cram" \) | head -n 1 || true
}
CRAM=$(resolve_cram_in_root "$SAMPLE_ID" "$CRAM_ROOT_1")
[[ -n "$CRAM" ]] || CRAM=$(resolve_cram_in_root "$SAMPLE_ID" "$CRAM_ROOT_2")
[[ -n "$CRAM" ]] || { echo "ERROR: CRAM not found for $SAMPLE_ID" >&2; exit 1; }

if [[ -f "${CRAM}.crai" ]]; then
  CRAI="${CRAM}.crai"
elif [[ -f "${CRAM%.cram}.cram.crai" ]]; then
  CRAI="${CRAM%.cram}.cram.crai"
else
  echo "ERROR: CRAI not found for $CRAM" >&2
  exit 1
fi

WGS_REF="${WHOLE_REF_DIR}/${REF_SPECIES}.fasta"
NUCLEAR_REF="${NUCLEAR_ONLY_REF_DIR}/${REF_SPECIES}.nuclear_only.fa"
[[ -s "$WGS_REF" ]] || { echo "ERROR: WGS ref missing: $WGS_REF" >&2; exit 1; }
[[ -s "$NUCLEAR_REF" ]] || { echo "ERROR: nuclear ref missing: $NUCLEAR_REF" >&2; exit 1; }

samtools faidx "$WGS_REF" >/dev/null 2>&1 || true
MT_CONTIG="${MT_CONTIG:-chrM}"
MT_LENGTH=$(awk -v mt="$MT_CONTIG" '$1==mt{print $2; exit}' "${WGS_REF}.fai")
[[ -n "$MT_LENGTH" ]] || { echo "ERROR: mt contig $MT_CONTIG not found in ${WGS_REF}.fai" >&2; exit 1; }

SAMPLE_DISCOVERY_OUTDIR="${DISCOVERY_OUTROOT}/${SAMPLE_ID}"

# Do not skip solely on existing discovery outputs:
# run_numt_end2end.sh handles skip logic using final highconf output status,
# so best-hit can still run when discovery BED/TSV already exist.

TMP_CFG=$(mktemp "${TMPDIR:-/tmp}/${SAMPLE_ID}.numtcfg.XXXXXX")
cp "$CONFIG" "$TMP_CFG"
cat >> "$TMP_CFG" <<CFG
SAMPLE=${SAMPLE_ID}
INPUT_BAM_CRAM=${CRAM}
INPUT_INDEX=${CRAI}
WGS_REF=${WGS_REF}
NUCLEAR_REF=${NUCLEAR_REF}
MT_CONTIG=${MT_CONTIG}
MT_LENGTH=${MT_LENGTH}
DISCOVERY_OUTDIR=${SAMPLE_DISCOVERY_OUTDIR}
SAMPLES_TSV=${SAMPLES_TSV}
WHOLE_REF_DIR=${WHOLE_REF_DIR}
CHRM_REF_DIR=${CHRM_REF_DIR}
BESTHIT_OUTDIR=${BESTHIT_OUTDIR}
CFG

bash "${REPO_DIR}/run_numt_end2end.sh" --config "$TMP_CFG"
rm -f "$TMP_CFG"
