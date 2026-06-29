#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}"
source "${REPO_DIR}/load_numt_modules.sh"

########################################
# Usage
########################################
usage() {
  cat <<EOF
Usage:
  bash run_numt_discovery.sh \\
    --input sample.cram \\
    --index sample.cram.crai \\
    --sample SAMPLE1 \\
    --mt-contig chrM \\
    --mt-length 16569 \\
    --wgs-ref ref.fa \\
    --nuclear-ref nuclear_only.fa \\
    --outdir results/SAMPLE1 \\
    [--threads 8] \\
    [--min-mapq 20] \\
    [--min-depth 3] \\
    [--min-reads 5] \\
    [--min-len 100] \\
    [--merge-gap 50] \\
    [--pad 500] \
    [--input-alt sample_alt.cram] \
    [--index-alt sample_alt.cram.crai]

Required:
  --input         Input BAM/CRAM
  --index         BAM/CRAM index
  --sample        Sample name
  --mt-contig     mt contig name, e.g. chrM
  --mt-length     mt length, e.g. 16569
  --wgs-ref       Full reference fasta (required for CRAM)
  --nuclear-ref   Nuclear-only fasta used for remapping
  --outdir        Output directory
  --input-alt     Optional alternate BAM/CRAM path if --input is missing
  --index-alt     Optional alternate BAM/CRAM index path if --index is missing

Notes:
  1. nuclear-only fasta should NOT contain chrM.
  2. bwa index should already be built for nuclear-only fasta.
  3. samtools faidx should already be built for all fasta files.
EOF
}

########################################
# Defaults
########################################
THREADS=8
MIN_MAPQ=20
MIN_DEPTH=3
MIN_READS=5
MIN_LEN=100
MERGE_GAP=50
PAD=500
FORCE_RERUN=0

########################################
# Parse args
########################################
INPUT=""
INDEX=""
SAMPLE=""
MT_CONTIG=""
MT_LENGTH=""
WGS_REF=""
NUCLEAR_REF=""
OUTDIR=""
CONFIG=""
INPUT_ALT=""
INDEX_ALT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) INPUT="$2"; shift 2 ;;
    --index) INDEX="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --mt-contig) MT_CONTIG="$2"; shift 2 ;;
    --mt-length) MT_LENGTH="$2"; shift 2 ;;
    --wgs-ref) WGS_REF="$2"; shift 2 ;;
    --nuclear-ref) NUCLEAR_REF="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --config) CONFIG="$2"; shift 2 ;;
    --input-alt) INPUT_ALT="$2"; shift 2 ;;
    --index-alt) INDEX_ALT="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
    --min-depth) MIN_DEPTH="$2"; shift 2 ;;
    --min-reads) MIN_READS="$2"; shift 2 ;;
    --min-len) MIN_LEN="$2"; shift 2 ;;
    --merge-gap) MERGE_GAP="$2"; shift 2 ;;
    --pad) PAD="$2"; shift 2 ;;
    --force-rerun) FORCE_RERUN=1; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done


if [[ -n "$CONFIG" ]]; then
  [[ -s "$CONFIG" ]] || { echo "ERROR: --config file not found: $CONFIG" >&2; exit 1; }
  # shellcheck disable=SC1090
  source "$CONFIG"

  INPUT="${INPUT:-${INPUT_BAM_CRAM:-$INPUT}}"
  INDEX="${INDEX:-${INPUT_INDEX:-$INDEX}}"
  SAMPLE="${SAMPLE:-$SAMPLE}"
  MT_CONTIG="${MT_CONTIG:-$MT_CONTIG}"
  MT_LENGTH="${MT_LENGTH:-$MT_LENGTH}"
  WGS_REF="${WGS_REF:-$WGS_REF}"
  NUCLEAR_REF="${NUCLEAR_REF:-$NUCLEAR_REF}"
  OUTDIR="${OUTDIR:-${DISCOVERY_OUTDIR:-$OUTDIR}}"
  INPUT_ALT="${INPUT_ALT:-${INPUT_BAM_CRAM_ALT:-$INPUT_ALT}}"
  INDEX_ALT="${INDEX_ALT:-${INPUT_INDEX_ALT:-$INDEX_ALT}}"

  THREADS="${THREADS:-$THREADS}"
  MIN_MAPQ="${MIN_MAPQ_DISCOVERY:-$MIN_MAPQ}"
  MIN_DEPTH="${MIN_DEPTH:-$MIN_DEPTH}"
  MIN_READS="${MIN_READS:-$MIN_READS}"
  MIN_LEN="${MIN_LEN:-$MIN_LEN}"
  MERGE_GAP="${MERGE_GAP_DISCOVERY:-$MERGE_GAP}"
  PAD="${PAD:-$PAD}"
fi

# If primary input/index missing, fallback to alternate paths if provided
if [[ -n "$INPUT" && ! -e "$INPUT" && -n "$INPUT_ALT" && -e "$INPUT_ALT" ]]; then
  echo "[INFO] Primary --input not found. Using alternate input path: $INPUT_ALT"
  INPUT="$INPUT_ALT"
fi

if [[ -n "$INDEX" && ! -e "$INDEX" && -n "$INDEX_ALT" && -e "$INDEX_ALT" ]]; then
  echo "[INFO] Primary --index not found. Using alternate index path: $INDEX_ALT"
  INDEX="$INDEX_ALT"
fi

########################################
# Check inputs
########################################
[[ -n "$INPUT" ]] || { echo "ERROR: --input required" >&2; exit 1; }
[[ -n "$INDEX" ]] || { echo "ERROR: --index required" >&2; exit 1; }
[[ -n "$SAMPLE" ]] || { echo "ERROR: --sample required" >&2; exit 1; }
[[ -n "$MT_CONTIG" ]] || { echo "ERROR: --mt-contig required" >&2; exit 1; }
[[ -n "$MT_LENGTH" ]] || { echo "ERROR: --mt-length required" >&2; exit 1; }
[[ -n "$WGS_REF" ]] || { echo "ERROR: --wgs-ref required" >&2; exit 1; }
[[ -n "$NUCLEAR_REF" ]] || { echo "ERROR: --nuclear-ref required" >&2; exit 1; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir required" >&2; exit 1; }
[[ -s "$INPUT" ]] || { echo "ERROR: input file not found or empty: $INPUT" >&2; exit 1; }
[[ -s "$INDEX" ]] || { echo "ERROR: index file not found or empty: $INDEX" >&2; exit 1; }

mkdir -p "$OUTDIR"/{tmp,logs,intermediate}

LOG="$OUTDIR/logs/${SAMPLE}.numt_discovery.log"
exec > >(tee -a "$LOG") 2>&1

echo "[$(date)] Starting NUMT discovery for $SAMPLE"
echo "Input: $INPUT"
echo "Index: $INDEX"
echo "mt contig: $MT_CONTIG"
echo "mt length: $MT_LENGTH"
echo "WGS ref: $WGS_REF"
echo "Nuclear ref: $NUCLEAR_REF"
echo "Outdir: $OUTDIR"

########################################
# Files
########################################
MT_BED="$OUTDIR/intermediate/${SAMPLE}.mt_region.bed"
MT_FETCH_BAM="$OUTDIR/intermediate/${SAMPLE}.mt_fetch_pairs.bam"
MT_FETCH_BAI="$OUTDIR/intermediate/${SAMPLE}.mt_fetch_pairs.bam.bai"

R1_FQ="$OUTDIR/intermediate/${SAMPLE}.R1.fastq.gz"
R2_FQ="$OUTDIR/intermediate/${SAMPLE}.R2.fastq.gz"
SINGLE_FQ="$OUTDIR/intermediate/${SAMPLE}.single.fastq.gz"

NUC_BAM="$OUTDIR/intermediate/${SAMPLE}.mt_related_to_nuclear.sorted.bam"
NUC_BAI="$OUTDIR/intermediate/${SAMPLE}.mt_related_to_nuclear.sorted.bam.bai"
NUC_BAM_TMP="$OUTDIR/intermediate/${SAMPLE}.mt_related_to_nuclear.sorted.bam.tmp"
NUC_BAI_TMP="$OUTDIR/intermediate/${SAMPLE}.mt_related_to_nuclear.sorted.bam.tmp.bai"

BED_OUT="$OUTDIR/${SAMPLE}.numt_candidates.bed"
TSV_OUT="$OUTDIR/${SAMPLE}.numt_candidates.tsv"
DONE_OUT="$OUTDIR/${SAMPLE}.numt_discovery.done"
HIGHCONF_OUT="$OUTDIR/${SAMPLE}.highconf_numt.bed"

########################################
# Resume/skip guard
########################################
if [[ "$FORCE_RERUN" -ne 1 && -f "$DONE_OUT" && -e "$BED_OUT" && -s "$TSV_OUT" ]]; then
  echo "[$(date)] Existing completion marker detected; sample appears completed."
  echo "[$(date)] Skipping discovery for $SAMPLE"
  echo "DONE: $DONE_OUT"
  echo "BED: $BED_OUT"
  echo "TSV: $TSV_OUT"
  exit 0
fi

########################################
# Step 1. create mt BED
########################################
echo -e "${MT_CONTIG}\t1\t${MT_LENGTH}" > "$MT_BED"

########################################
# Step 2. fetch chrM reads + mates from WGS
########################################
# -P / --fetch-pairs fetches mates for reads overlapping region
# For CRAM, -T ref.fa is required
echo "[$(date)] Step 2: extracting mt reads + mates"

samtools view \
  -@ "$THREADS" \
  -T "$WGS_REF" \
  -P \
  -L "$MT_BED" \
  -b \
  -h \
  "$INPUT" > "$MT_FETCH_BAM"

samtools index -@ "$THREADS" "$MT_FETCH_BAM" "$MT_FETCH_BAI"

########################################
# Step 3. BAM -> FASTQ
########################################
echo "[$(date)] Step 3: converting BAM to FASTQ"

samtools fastq \
  -@ "$THREADS" \
  -1 "$R1_FQ" \
  -2 "$R2_FQ" \
  -0 /dev/null \
  -s "$SINGLE_FQ" \
  -n \
  "$MT_FETCH_BAM"

########################################
# Step 4. align to nuclear-only reference
########################################
echo "[$(date)] Step 4: remapping to nuclear-only reference"

# avoid stale/truncated outputs from interrupted previous runs
PAIR_BAM_TMP="$OUTDIR/intermediate/${SAMPLE}.mt_related_to_nuclear.paired.sorted.bam"
SINGLE_BAM_TMP="$OUTDIR/intermediate/${SAMPLE}.mt_related_to_nuclear.single.sorted.bam"
rm -f "$NUC_BAM" "$NUC_BAI" "$NUC_BAM_TMP" "$NUC_BAI_TMP" "$PAIR_BAM_TMP" "$SINGLE_BAM_TMP"

HAS_PAIRED=0
HAS_SINGLE=0
[[ -s "$R1_FQ" && -s "$R2_FQ" ]] && HAS_PAIRED=1
[[ -s "$SINGLE_FQ" ]] && HAS_SINGLE=1

if [[ "$HAS_PAIRED" -eq 0 && "$HAS_SINGLE" -eq 0 ]]; then
  echo "[$(date)] WARNING: no reads available for remapping (paired and singleton FASTQ are empty)."
  : > "$BED_OUT"
  : > "$HIGHCONF_OUT"
  echo -e "sample\tchrom\tstart0\tend0\tlength\tnreads\tmean_mapq\tpadded_start0\tpadded_end0\tpadded_length" > "$TSV_OUT"
  touch "$DONE_OUT"
  echo "[$(date)] Done (no candidates)."
  exit 0
fi

PAIR_SAM_TMP="$OUTDIR/intermediate/${SAMPLE}.mt_related_to_nuclear.paired.sam"
SINGLE_SAM_TMP="$OUTDIR/intermediate/${SAMPLE}.mt_related_to_nuclear.single.sam"
rm -f "$PAIR_SAM_TMP" "$SINGLE_SAM_TMP"

if [[ "$HAS_PAIRED" -eq 1 ]]; then
  bwa mem \
    -t "$THREADS" \
    -K 100000000 \
    "$NUCLEAR_REF" \
    "$R1_FQ" "$R2_FQ" > "$PAIR_SAM_TMP"
  if [[ -s "$PAIR_SAM_TMP" ]] && samtools view -H "$PAIR_SAM_TMP" >/dev/null 2>&1; then
    samtools sort -@ "$THREADS" -o "$PAIR_BAM_TMP" "$PAIR_SAM_TMP"
  else
    echo "[$(date)] WARNING: paired remap produced no valid SAM output; skipping paired BAM."
    HAS_PAIRED=0
  fi
fi

if [[ "$HAS_SINGLE" -eq 1 ]]; then
  bwa mem \
    -t "$THREADS" \
    -K 100000000 \
    "$NUCLEAR_REF" \
    "$SINGLE_FQ" > "$SINGLE_SAM_TMP"
  if [[ -s "$SINGLE_SAM_TMP" ]] && samtools view -H "$SINGLE_SAM_TMP" >/dev/null 2>&1; then
    samtools sort -@ "$THREADS" -o "$SINGLE_BAM_TMP" "$SINGLE_SAM_TMP"
  else
    echo "[$(date)] WARNING: singleton remap produced no valid SAM output; skipping singleton BAM."
    HAS_SINGLE=0
  fi
fi

rm -f "$PAIR_SAM_TMP" "$SINGLE_SAM_TMP"

if [[ "$HAS_PAIRED" -eq 0 && "$HAS_SINGLE" -eq 0 ]]; then
  echo "[$(date)] WARNING: no valid remap outputs were produced."
  : > "$BED_OUT"
  : > "$HIGHCONF_OUT"
  echo -e "sample\tchrom\tstart0\tend0\tlength\tnreads\tmean_mapq\tpadded_start0\tpadded_end0\tpadded_length" > "$TSV_OUT"
  touch "$DONE_OUT"
  echo "[$(date)] Done (no candidates)."
  exit 0
fi

if [[ "$HAS_PAIRED" -eq 1 && "$HAS_SINGLE" -eq 1 ]]; then
  samtools merge -@ "$THREADS" -f "$NUC_BAM_TMP" "$PAIR_BAM_TMP" "$SINGLE_BAM_TMP"
elif [[ "$HAS_PAIRED" -eq 1 ]]; then
  mv "$PAIR_BAM_TMP" "$NUC_BAM_TMP"
else
  mv "$SINGLE_BAM_TMP" "$NUC_BAM_TMP"
fi

# validate BAM integrity before moving to final path
samtools quickcheck -v "$NUC_BAM_TMP"

NUC_REC_COUNT=$(samtools view -c "$NUC_BAM_TMP")
if [[ "$NUC_REC_COUNT" -eq 0 ]]; then
  echo "[$(date)] WARNING: remapped BAM has 0 records."
  echo "[$(date)] No mt-related nuclear alignments found for ${SAMPLE}."

  mv "$NUC_BAM_TMP" "$NUC_BAM"
   : > "$BED_OUT"
  : > "$HIGHCONF_OUT"
  echo -e "sample	chrom	start0	end0	length	nreads	mean_mapq	padded_start0	padded_end0	padded_length" > "$TSV_OUT"
  touch "$DONE_OUT"

  echo "[$(date)] Done (no candidates)."
  echo "BED: $BED_OUT (empty)"
  echo "HIGHCONF BED: $HIGHCONF_OUT (empty)"
  echo "TSV: $TSV_OUT (header only)"
  exit 0
fi

samtools index -@ "$THREADS" "$NUC_BAM_TMP" "$NUC_BAI_TMP"
samtools quickcheck -v "$NUC_BAM_TMP"

mv "$NUC_BAM_TMP" "$NUC_BAM"
mv "$NUC_BAI_TMP" "$NUC_BAI"

echo "[$(date)] Step 4 check: BAM + BAI ready"

########################################
# Step 5. discover sink loci
########################################
echo "[$(date)] Step 5: discovering candidate NUMT sink loci"

python3 discover_numt_sinks.py \
  --bam "$NUC_BAM" \
  --sample "$SAMPLE" \
  --min-mapq "$MIN_MAPQ" \
  --min-depth "$MIN_DEPTH" \
  --min-reads "$MIN_READS" \
  --min-len "$MIN_LEN" \
  --merge-gap "$MERGE_GAP" \
  --pad "$PAD" \
  --bed-out "$BED_OUT" \
  --tsv-out "$TSV_OUT"

touch "$DONE_OUT"
echo "[$(date)] Done."
echo "DONE: $DONE_OUT"
echo "BED: $BED_OUT"
echo "TSV: $TSV_OUT"
