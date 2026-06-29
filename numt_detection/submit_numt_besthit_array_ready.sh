#!/usr/bin/env bash
#SBATCH --job-name=numt_besthit
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/numt_besthit_%A_%a.log

set -euo pipefail

INPUT_DIR="/nfs/roberts/pi/pi_njl27/lt692/numt_discovery/numt_candidates_bed_strict"
SAMPLES_TSV="/nfs/roberts/pi/pi_njl27/lt692/numt_discovery/human_sample.txt"
WHOLE_REF_DIR="/nfs/roberts/pi/pi_njl27/lt692/primate_ref_files/whole_genomes_final"
CHRM_REF_DIR="/nfs/roberts/pi/pi_njl27/lt692/primate_ref_files/Refseq_chrM"
OUTDIR="/nfs/roberts/pi/pi_njl27/lt692/numt_discovery/primate_results_numt_besthit_strict"
CONCURRENT=20
MERGE_GAP=50
MIN_PIDENT=90
MIN_ALIGN_LENGTH=50
MIN_MAPQ=20
MIN_NREADS=5
RECURSIVE_FIND=1
WORKER_SCRIPT="/nfs/roberts/pi/pi_njl27/lt692/numt_discovery/process_numt_candidates_one_ready.sh"

usage() {
  cat <<USAGE
Usage:
  $(basename "$0") [optional overrides]

Overrides:
  --input-dir         Default: ${INPUT_DIR}
  --samples-tsv       Default: ${SAMPLES_TSV}
  --whole-ref-dir     Default: ${WHOLE_REF_DIR}
  --chrm-ref-dir      Default: ${CHRM_REF_DIR}
  --outdir            Default: ${OUTDIR}
  --concurrent        Default: ${CONCURRENT}
  --merge-gap         Default: ${MERGE_GAP}
  --min-pident        Default: ${MIN_PIDENT}
  --min-align-length  Default: ${MIN_ALIGN_LENGTH}
  --min-mapq          Default: ${MIN_MAPQ}
  --min-nreads        Default: ${MIN_NREADS}
  --no-recursive      Only scan INPUT_DIR top level
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir) INPUT_DIR="$2"; shift 2 ;;
    --samples-tsv) SAMPLES_TSV="$2"; shift 2 ;;
    --whole-ref-dir) WHOLE_REF_DIR="$2"; shift 2 ;;
    --chrm-ref-dir) CHRM_REF_DIR="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --concurrent) CONCURRENT="$2"; shift 2 ;;
    --merge-gap) MERGE_GAP="$2"; shift 2 ;;
    --min-pident) MIN_PIDENT="$2"; shift 2 ;;
    --min-align-length) MIN_ALIGN_LENGTH="$2"; shift 2 ;;
    --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
    --min-nreads) MIN_NREADS="$2"; shift 2 ;;
    --no-recursive) RECURSIVE_FIND=0; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -d "$INPUT_DIR" ]] || { echo "ERROR: input dir not found: $INPUT_DIR" >&2; exit 1; }
[[ -s "$SAMPLES_TSV" ]] || { echo "ERROR: samples TSV not found: $SAMPLES_TSV" >&2; exit 1; }
[[ -d "$WHOLE_REF_DIR" ]] || { echo "ERROR: whole-ref dir not found: $WHOLE_REF_DIR" >&2; exit 1; }
[[ -d "$CHRM_REF_DIR" ]] || { echo "ERROR: chrM-ref dir not found: $CHRM_REF_DIR" >&2; exit 1; }
[[ -x "$WORKER_SCRIPT" ]] || { echo "ERROR: worker script not executable: $WORKER_SCRIPT" >&2; exit 1; }

mkdir -p "$OUTDIR" logs

# 按 samples.tsv 行数决定 array 大小
N=$(wc -l < "$SAMPLES_TSV")
[[ "$N" -gt 0 ]] || { echo "ERROR: samples TSV is empty: $SAMPLES_TSV" >&2; exit 1; }

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "Found $N samples in $SAMPLES_TSV"
  echo "Input dir : $INPUT_DIR"
  echo "Output dir: $OUTDIR"
  echo "Submitting array 1-$N with %$CONCURRENT concurrency"

  sbatch --array=1-${N}%${CONCURRENT} "$0" \
    --input-dir "$INPUT_DIR" \
    --samples-tsv "$SAMPLES_TSV" \
    --whole-ref-dir "$WHOLE_REF_DIR" \
    --chrm-ref-dir "$CHRM_REF_DIR" \
    --outdir "$OUTDIR" \
    --concurrent "$CONCURRENT" \
    --merge-gap "$MERGE_GAP" \
    --min-pident "$MIN_PIDENT" \
    --min-align-length "$MIN_ALIGN_LENGTH" \
    --min-mapq "$MIN_MAPQ" \
    --min-nreads "$MIN_NREADS" \
    $([[ "$RECURSIVE_FIND" -eq 0 ]] && echo --no-recursive)
  exit 0
fi

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_TSV" || true)
[[ -n "$LINE" ]] || { echo "ERROR: no line found for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }

SAMPLE_ID=$(printf '%s\n' "$LINE" | cut -f1)
[[ -n "$SAMPLE_ID" ]] || { echo "ERROR: empty sample_id at task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }

# 根据 sample_id 找对应的 BED
if [[ "$RECURSIVE_FIND" -eq 1 ]]; then
  BED=$(find "$INPUT_DIR" -type f -name "${SAMPLE_ID}.numt_candidates.bed" | head -n 1 || true)
else
  BED=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "${SAMPLE_ID}.numt_candidates.bed" | head -n 1 || true)
fi

# 如果找不到 BED，跳过
if [[ -z "$BED" ]]; then
  echo "[$(date)] Skip ${SAMPLE_ID}: BED not found under $INPUT_DIR"
  exit 0
fi

# 如果 BED 为空，跳过
if [[ ! -s "$BED" ]]; then
  echo "[$(date)] Skip ${SAMPLE_ID}: empty BED -> $BED"
  exit 0
fi

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID} processing sample ${SAMPLE_ID}"
echo "[$(date)] BED: $BED"

bash "$WORKER_SCRIPT" \
  --bed "$BED" \
  --samples-tsv "$SAMPLES_TSV" \
  --whole-ref-dir "$WHOLE_REF_DIR" \
  --chrm-ref-dir "$CHRM_REF_DIR" \
  --outdir "$OUTDIR" \
  --merge-gap "$MERGE_GAP" \
  --min-pident "$MIN_PIDENT" \
  --min-align-length "$MIN_ALIGN_LENGTH" \
  --min-mapq "$MIN_MAPQ" \
  --min-nreads "$MIN_NREADS"