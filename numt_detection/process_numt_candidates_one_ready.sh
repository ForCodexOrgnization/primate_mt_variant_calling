#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}"
source "${REPO_DIR}/load_numt_modules.sh"

# ============================================================
# Ready-to-use defaults for your Yale HPC layout
# Can still be overridden by command-line arguments.
# ============================================================
DEFAULT_SAMPLES_TSV="/nfs/roberts/project/pi_njl27/lt692/primate_mito_calling/samples.tsv"
DEFAULT_WHOLE_REF_DIR="/nfs/roberts/pi/pi_njl27/lt692/primate_ref_files/whole_genomes_final"
DEFAULT_CHRM_REF_DIR="/nfs/roberts/pi/pi_njl27/lt692/primate_ref_files/Refseq_chrM"
DEFAULT_OUTDIR="/nfs/roberts/pi/pi_njl27/lt692/primate_results_numt_besthit"
DEFAULT_MERGE_GAP=50
DEFAULT_MIN_PIDENT=90
DEFAULT_MIN_ALIGN_LENGTH=50
DEFAULT_MIN_MAPQ=30
DEFAULT_MIN_NREADS=5
# ============================================================

usage() {
  cat <<USAGE
Usage:
  $(basename "$0") --bed SAMPLE.numt_candidates.bed [options]

Options:
  --bed               Input bed file, e.g. ERS12091914.numt_candidates.bed
  --samples-tsv       Default: ${DEFAULT_SAMPLES_TSV}
  --whole-ref-dir     Default: ${DEFAULT_WHOLE_REF_DIR}
  --chrm-ref-dir      Default: ${DEFAULT_CHRM_REF_DIR}
  --outdir            Default: ${DEFAULT_OUTDIR}
  --merge-gap         Default: ${DEFAULT_MERGE_GAP}
  --min-pident        Default: ${DEFAULT_MIN_PIDENT}
  --min-align-length  Default: ${DEFAULT_MIN_ALIGN_LENGTH}
  --min-mapq          Default: ${DEFAULT_MIN_MAPQ}
  --min-nreads        Default: ${DEFAULT_MIN_NREADS}
USAGE
}

BED=""
SAMPLES_TSV="$DEFAULT_SAMPLES_TSV"
WHOLE_REF_DIR="$DEFAULT_WHOLE_REF_DIR"
CHRM_REF_DIR="$DEFAULT_CHRM_REF_DIR"
OUTDIR="$DEFAULT_OUTDIR"
MERGE_GAP="$DEFAULT_MERGE_GAP"
MIN_PIDENT="$DEFAULT_MIN_PIDENT"
MIN_ALIGN_LENGTH="$DEFAULT_MIN_ALIGN_LENGTH"
MIN_MAPQ="$DEFAULT_MIN_MAPQ"
MIN_NREADS="$DEFAULT_MIN_NREADS"
CONFIG=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bed) BED="$2"; shift 2 ;;
    --config) CONFIG="$2"; shift 2 ;;
    --samples-tsv) SAMPLES_TSV="$2"; shift 2 ;;
    --whole-ref-dir) WHOLE_REF_DIR="$2"; shift 2 ;;
    --chrm-ref-dir) CHRM_REF_DIR="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --merge-gap) MERGE_GAP="$2"; shift 2 ;;
    --min-pident) MIN_PIDENT="$2"; shift 2 ;;
    --min-align-length) MIN_ALIGN_LENGTH="$2"; shift 2 ;;
    --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
    --min-nreads) MIN_NREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -n "$CONFIG" ]]; then
  [[ -s "$CONFIG" ]] || { echo "ERROR: --config file not found: $CONFIG" >&2; exit 1; }
  # shellcheck disable=SC1090
  source "$CONFIG"

  SAMPLES_TSV="${SAMPLES_TSV:-$SAMPLES_TSV}"
  WHOLE_REF_DIR="${WHOLE_REF_DIR:-$WHOLE_REF_DIR}"
  CHRM_REF_DIR="${CHRM_REF_DIR:-$CHRM_REF_DIR}"
  OUTDIR="${BESTHIT_OUTDIR:-$OUTDIR}"
  MERGE_GAP="${MERGE_GAP_BESTHIT:-$MERGE_GAP}"
  MIN_PIDENT="${MIN_PIDENT:-$MIN_PIDENT}"
  MIN_ALIGN_LENGTH="${MIN_ALIGN_LENGTH:-$MIN_ALIGN_LENGTH}"
  MIN_MAPQ="${MIN_MAPQ_BESTHIT:-$MIN_MAPQ}"
  MIN_NREADS="${MIN_NREADS:-$MIN_NREADS}"
fi

[[ -n "$BED" ]] || { echo "ERROR: missing --bed" >&2; usage; exit 1; }
[[ -e "$BED" ]] || { echo "ERROR: input BED not found: $BED" >&2; exit 1; }
[[ -s "$SAMPLES_TSV" ]] || { echo "ERROR: samples TSV not found: $SAMPLES_TSV" >&2; exit 1; }
[[ -d "$WHOLE_REF_DIR" ]] || { echo "ERROR: whole-ref dir not found: $WHOLE_REF_DIR" >&2; exit 1; }
[[ -d "$CHRM_REF_DIR" ]] || { echo "ERROR: chrM-ref dir not found: $CHRM_REF_DIR" >&2; exit 1; }

mkdir -p "$OUTDIR"

for exe in awk sort bedtools blastn makeblastdb; do
  command -v "$exe" >/dev/null 2>&1 || {
    echo "ERROR: required executable not found: $exe" >&2
    exit 1
  }
done

sample_id=$(basename "$BED")
sample_id=${sample_id%.numt_candidates.bed}

# sample_id in col1, ref_name in col3; support tab/comma/space-delimited sample sheets.
ref_name=$(awk -v s="$sample_id" '
function trim(x){ gsub(/^[[:space:]]+|[[:space:]]+$/, "", x); return x }
{
  line=$0
  gsub(/\r/, "", line)
  if (line ~ /^[[:space:]]*$/ || line ~ /^[[:space:]]*#/) next

  n=split(line, a, /[\t,]/)
  k=0
  for (i=1; i<=n; i++) {
    v=trim(a[i])
    if (v != "") {
      k++
      f[k]=v
      if (k==3) break
    }
  }

  if (k < 3) {
    n2=split(line, b, /[[:space:]]+/)
    k=0
    for (j=1; j<=n2; j++) {
      v2=trim(b[j])
      if (v2 != "") {
        k++
        f[k]=v2
        if (k==3) break
      }
    }
  }

  if (k >= 3 && f[1] == s) {
    print f[3]
    exit
  }
}
' "$SAMPLES_TSV")
[[ -n "$ref_name" ]] || {
  echo "ERROR: sample $sample_id not found in $SAMPLES_TSV" >&2
  exit 1
}

whole_ref="${WHOLE_REF_DIR}/${ref_name}.fasta"
chrm_ref="${CHRM_REF_DIR}/${ref_name}.fasta"

[[ -s "$whole_ref" ]] || { echo "ERROR: whole-genome fasta not found: $whole_ref" >&2; exit 1; }
[[ -s "$chrm_ref"  ]] || { echo "ERROR: chrM fasta not found: $chrm_ref" >&2; exit 1; }

prefix="${OUTDIR}/${sample_id}"
besthit_tsv="${prefix}.numt_vs_chrM.besthit.tsv"
highconf_bed="${prefix}.highconf_numt.bed"

if [[ ! -s "$BED" ]]; then
  echo "[INFO] Input BED is empty for ${sample_id}; writing empty highconf output and exiting."
  : > "$highconf_bed"
  echo -e "locus_id	nuclear_chrom	nuclear_start	nuclear_end	subject	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	plus_nreads	plus_meanMAPQ	has_numt_signal" > "$besthit_tsv"
  echo "[INFO] Wrote: $highconf_bed (empty)"
  echo "[INFO] Wrote: $besthit_tsv (header only)"
  exit 0
fi

tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/${sample_id}.numt_besthit.XXXXXX")
cleanup() {
  rm -rf "$tmpdir"
}
trap cleanup EXIT

summary_bed="${tmpdir}/${sample_id}.merged.summary.bed"
summary_tsv="${tmpdir}/${sample_id}.merged.summary.tsv"
named_bed="${tmpdir}/${sample_id}.merged.named.bed"
regions_fa="${tmpdir}/${sample_id}.merged.regions.fa"
blast_db_prefix="${tmpdir}/${sample_id}.chrMdb"
blast_raw_tsv="${tmpdir}/${sample_id}.blast.raw.tsv"
blast_best_tsv="${tmpdir}/${sample_id}.blast.best.tsv"

# ------------------------------------------------------------
# 1) Merge nearby candidate blocks and summarize per merged locus
#    Supports both:
#      BED3: chr start end
#      BED4+: chr start end info_with_nreads_meanMAPQ
# ------------------------------------------------------------

N_COL=$(awk '
  BEGIN{FS="[ \t]+"}
  $0 !~ /^#/ && NF > 0 {
    print NF
    exit
  }
' "$BED")

if [[ -z "${N_COL}" ]]; then
  echo "ERROR: cannot determine BED column number: $BED" >&2
  exit 1
fi

if [[ "$N_COL" -ge 4 ]]; then
  echo "[INFO] Input BED has >=4 columns. Will parse nreads / meanMAPQ from column 4."

  sort -k1,1 -k2,2n -k3,3n "$BED" \
    | bedtools merge -d "$MERGE_GAP" -c 4 -o collapse -i - \
    | awk -F'\t' -v OFS='\t' '
        function trim(x){ gsub(/^ +| +$/, "", x); return x }
        {
          chr=$1; start=$2; end=$3; collapsed=$4;
          n=split(collapsed, a, ",");

          sum_nreads=0;
          max_nreads=0;
          sum_mapq=0;
          n_mapq=0;

          for(i=1;i<=n;i++){
            block=trim(a[i]);
            nreads=0;
            mapq=0;

            m=split(block, parts, ";");
            for(j=1;j<=m;j++){
              if(parts[j] ~ /^nreads=/){
                split(parts[j], x, "=");
                nreads=x[2]+0;
              } else if(parts[j] ~ /^meanMAPQ=/){
                split(parts[j], y, "=");
                mapq=y[2]+0;
              }
            }

            sum_nreads += nreads;
            if(nreads > max_nreads) max_nreads=nreads;
            sum_mapq += mapq;
            n_mapq++;
          }

          mean_nreads = (n > 0 ? sum_nreads / n : 0);
          mean_mapq   = (n_mapq > 0 ? sum_mapq / n_mapq : 0);
          locus_id    = chr ":" start "-" end;

          print chr, start, end, n, sum_nreads, max_nreads, mean_nreads, mean_mapq, locus_id;
        }
      ' > "$summary_bed"

else
  echo "[INFO] Input BED has 3 columns. Treating as reference-level NUMT BED."
  echo "[INFO] nreads / meanMAPQ will be set to NA and will not be used for filtering."

  sort -k1,1 -k2,2n -k3,3n "$BED" \
    | bedtools merge -d "$MERGE_GAP" -i - \
    | awk -F'\t' -v OFS='\t' '
        {
          chr=$1; start=$2; end=$3;
          locus_id = chr ":" start "-" end;

          # columns:
          # chr start end n_blocks sum_nreads max_nreads mean_nreads mean_MAPQ locus_id
          print chr, start, end, 1, "NA", "NA", "NA", "NA", locus_id;
        }
      ' > "$summary_bed"
fi

if [[ ! -s "$summary_bed" ]]; then
  echo "ERROR: summary_bed was not generated or is empty: $summary_bed" >&2
  exit 1
fi

# Always generate summary_tsv from summary_bed
awk -F'\t' -v OFS='\t' '
  BEGIN{
    print "nuclear_chrom","nuclear_start","nuclear_end","n_blocks","sum_nreads","max_nreads","mean_nreads","mean_MAPQ","locus_id"
  }
  {
    print $1,$2,$3,$4,$5,$6,$7,$8,$9
  }
' "$summary_bed" > "$summary_tsv"

if [[ ! -s "$summary_tsv" ]]; then
  echo "ERROR: summary_tsv was not generated or is empty: $summary_tsv" >&2
  exit 1
fi

# ------------------------------------------------------------
# 2) IMPORTANT FIX:
#    getfasta -nameOnly uses BED column 4 as sequence name.
#    So create a BED4 with locus_id in column 4.
# ------------------------------------------------------------
awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$9}' "$summary_bed" > "$named_bed"

bedtools getfasta \
  -fi "$whole_ref" \
  -bed "$named_bed" \
  -nameOnly \
  -fo "$regions_fa"

# ------------------------------------------------------------
# 3) BLAST merged nuclear loci against species chrM reference
# ------------------------------------------------------------
makeblastdb -in "$chrm_ref" -dbtype nucl -out "$blast_db_prefix" >/dev/null

blastn \
  -query "$regions_fa" \
  -db "$blast_db_prefix" \
  -task blastn \
  -dust no \
  -evalue 1e-5 \
  -max_target_seqs 5 \
  -max_hsps 1 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
  > "$blast_raw_tsv"

# Keep the best hit per locus
awk -F'\t' -v OFS='\t' '
  {
    q    = $1;
    bits = $12 + 0;
    len  = $4  + 0;

    if(!(q in seen) || bits > best_bits[q] || (bits == best_bits[q] && len > best_len[q])){
      seen[q] = 1;
      best_bits[q] = bits;
      best_len[q] = len;
      line[q] = $0;
    }
  }
  END{
    for(q in line) print line[q]
  }
' "$blast_raw_tsv" | sort -k1,1 > "$blast_best_tsv"

# ------------------------------------------------------------
# 4) Join BLAST best hit back to summary table by locus_id
# ------------------------------------------------------------
awk -F'\t' -v OFS='\t' '
  NR==FNR {
    hit[$1] = $0;
    next;
  }
  FNR==1 {
    print $0, "chrM_hit","pident","align_length","mismatch","gapopen","qstart","qend","chrM_start","chrM_end","evalue","bitscore";
    next;
  }
  {
    locus = $9;
    if(locus in hit){
      split(hit[locus], a, FS);
      print $0, a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12];
    } else {
      print $0, ".", ".", ".", ".", ".", ".", ".", ".", ".", ".", ".";
    }
  }
' "$blast_best_tsv" "$summary_tsv" > "$besthit_tsv"

# ------------------------------------------------------------
# 5) Filter to high-confidence NUMTs
# ------------------------------------------------------------
awk -F'\t' -v OFS='\t' \
    -v min_pid="$MIN_PIDENT" \
    -v min_alen="$MIN_ALIGN_LENGTH" \
    -v min_mapq="$MIN_MAPQ" \
    -v min_nreads="$MIN_NREADS" '
NR==1{
  for(i=1;i<=NF;i++) c[$i]=i
  next
}
{
  chr   = $(c["nuclear_chrom"])
  start = $(c["nuclear_start"])
  end   = $(c["nuclear_end"])
  hit   = $(c["chrM_hit"])
  pid   = $(c["pident"])
  alen  = $(c["align_length"])
  mapq  = $(c["mean_MAPQ"])
  nread = $(c["sum_nreads"])

  if(hit=="." || hit=="") next
  if(pid+0   < min_pid) next
  if(alen+0  < min_alen) next
  if(mapq!="NA" && mapq!="." && mapq!="" && mapq+0 < min_mapq) next
  if(nread!="NA" && nread!="." && nread!="" && nread+0 < min_nreads) next
  if(end+0 <= start+0) next

  print chr, start, end
}
' "$besthit_tsv" \
  | sort -k1,1 -k2,2n -k3,3n \
  | bedtools merge -i - > "$highconf_bed"

echo "[DONE] $sample_id"
echo "  ref_name : $ref_name"
echo "  besthit  : $besthit_tsv"
echo "  highconf : $highconf_bed"
