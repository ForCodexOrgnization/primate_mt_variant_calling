#!/usr/bin/env bash
set -euo pipefail

########################################
# Usage
########################################
usage() {
  cat <<EOF
Usage:
  $(basename "$0") <numt_candidates.bed> <reference.fa> <out_prefix>

Example:
  $(basename "$0") ERS12091914.numt_candidates.bed /path/to/reference.fa ERS12091914

Input BED format:
  chrom  start  end  sample;nreads=25;meanMAPQ=31.12

Final output:
  <out_prefix>.numt_vs_chrM.besthit.tsv

Intermediate files:
  <out_prefix>.numt_candidates.sorted.bed
  <out_prefix>.numt_candidates.merged.bed
  <out_prefix>.numt_candidates.merged.summary.tsv
  <out_prefix>.numt_candidates.merged.fa
  <out_prefix>.chrM.fa
  <out_prefix>.numt_vs_chrM.tsv
  <out_prefix>.numt_vs_chrM.filtered.tsv
  <out_prefix>.numt_vs_chrM.besthit.raw.tsv
EOF
  exit 1
}

########################################
# Args
########################################
if [[ $# -ne 3 ]]; then
  usage
fi

BED_IN="$1"
REF_FA="$2"
OUT_PREFIX="$3"

########################################
# Parameters
########################################
MT_CONTIG="chrM"
MIN_ALIGN_LEN=100
MIN_PIDENT=80

SORTED_BED="${OUT_PREFIX}.numt_candidates.sorted.bed"
MERGED_BED="${OUT_PREFIX}.numt_candidates.merged.bed"
MERGED_SUMMARY="${OUT_PREFIX}.numt_candidates.merged.summary.tsv"
MERGED_FA="${OUT_PREFIX}.numt_candidates.merged.fa"
CHRM_FA="${OUT_PREFIX}.chrM.fa"
BLAST_OUT="${OUT_PREFIX}.numt_vs_chrM.tsv"
FILTERED_OUT="${OUT_PREFIX}.numt_vs_chrM.filtered.tsv"
BESTHIT_RAW="${OUT_PREFIX}.numt_vs_chrM.besthit.raw.tsv"
BESTHIT_FINAL="${OUT_PREFIX}.numt_vs_chrM.besthit.tsv"
TMP_BODY="${OUT_PREFIX}.numt_vs_chrM.besthit.body.tsv"



module load SAMtools/1.21-GCC-13.3.0
module load BEDTools/2.31.1-GCC-13.3.0
module load miniconda/24.11.3

# source ~/miniconda3/etc/profile.d/conda.sh
conda activate blast_env


########################################
# Check dependencies
########################################
for cmd in samtools bedtools makeblastdb blastn awk sort cut grep; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "[ERROR] Required command not found: $cmd" >&2
    exit 1
  fi
done

########################################
# Check input files
########################################
if [[ ! -s "$BED_IN" ]]; then
  echo "[ERROR] BED input not found or empty: $BED_IN" >&2
  exit 1
fi

if [[ ! -s "$REF_FA" ]]; then
  echo "[ERROR] Reference fasta not found or empty: $REF_FA" >&2
  exit 1
fi

########################################
# Index reference if needed
########################################
if [[ ! -s "${REF_FA}.fai" ]]; then
  echo "[INFO] Indexing reference fasta..."
  samtools faidx "$REF_FA"
fi

########################################
# Check chrM exists
########################################
if ! cut -f1 "${REF_FA}.fai" | grep -Fxq "$MT_CONTIG"; then
  echo "[ERROR] mt contig '$MT_CONTIG' not found in reference fasta index: ${REF_FA}.fai" >&2
  exit 1
fi

########################################
# 1. Sort original BED
########################################
echo "[INFO] Sorting BED..."
sort -k1,1 -k2,2n "$BED_IN" > "$SORTED_BED"

########################################
# 2. Merge overlapping intervals
########################################
echo "[INFO] Merging overlapping intervals..."
cut -f1-3 "$SORTED_BED" | bedtools merge -i - > "$MERGED_BED"

if [[ ! -s "$MERGED_BED" ]]; then
  echo "[ERROR] Merged BED is empty." >&2
  exit 1
fi

########################################
# 3. Summarize merged loci
########################################
echo "[INFO] Summarizing merged loci..."
awk '
BEGIN{
  OFS="\t";
}
function parse_nreads(info,    a,n,i){
  n=split(info,a,";");
  for(i=1;i<=n;i++){
    if(a[i] ~ /^nreads=/){
      sub(/^nreads=/,"",a[i]);
      return a[i]+0;
    }
  }
  return 0;
}
function parse_mapq(info,    a,n,i){
  n=split(info,a,";");
  for(i=1;i<=n;i++){
    if(a[i] ~ /^meanMAPQ=/){
      sub(/^meanMAPQ=/,"",a[i]);
      return a[i]+0;
    }
  }
  return 0;
}
FNR==NR{
  merged_chr[++m]=$1;
  merged_start[m]=$2;
  merged_end[m]=$3;
  next;
}
{
  chr=$1; start=$2; end=$3; info=$4;
  nreads=parse_nreads(info);
  mapq=parse_mapq(info);

  for(i=1;i<=m;i++){
    # overlap with merged interval
    if(chr==merged_chr[i] && start < merged_end[i] && end > merged_start[i]){
      key=merged_chr[i] FS merged_start[i] FS merged_end[i];
      n_blocks[key]++;
      sum_nreads[key]+=nreads;
      sum_mapq[key]+=mapq;
      if(!(key in max_nreads) || nreads > max_nreads[key]) max_nreads[key]=nreads;
      break;
    }
  }
}
END{
  print "nuclear_chrom","nuclear_start","nuclear_end","n_blocks","sum_nreads","max_nreads","mean_nreads","mean_MAPQ","locus_id";
  for(i=1;i<=m;i++){
    key=merged_chr[i] FS merged_start[i] FS merged_end[i];
    nb=(key in n_blocks ? n_blocks[key] : 0);
    sn=(key in sum_nreads ? sum_nreads[key] : 0);
    mx=(key in max_nreads ? max_nreads[key] : 0);
    mn=(nb>0 ? sn/nb : 0);
    mq=(nb>0 ? sum_mapq[key]/nb : 0);
    locus=merged_chr[i] ":" merged_start[i] "-" merged_end[i];
    print merged_chr[i], merged_start[i], merged_end[i], nb, sn, mx, mn, mq, locus;
  }
}
' "$MERGED_BED" "$SORTED_BED" > "$MERGED_SUMMARY"

########################################
# 4. Extract merged nuclear loci fasta
########################################
echo "[INFO] Extracting merged nuclear loci fasta..."
bedtools getfasta \
  -fi "$REF_FA" \
  -bed "$MERGED_BED" \
  -fo "$MERGED_FA"

if [[ ! -s "$MERGED_FA" ]]; then
  echo "[ERROR] Failed to generate merged fasta: $MERGED_FA" >&2
  exit 1
fi

########################################
# 5. Extract chrM fasta from same reference
########################################
echo "[INFO] Extracting ${MT_CONTIG} from reference..."
samtools faidx "$REF_FA" "$MT_CONTIG" > "$CHRM_FA"

if [[ ! -s "$CHRM_FA" ]]; then
  echo "[ERROR] Failed to generate chrM fasta: $CHRM_FA" >&2
  exit 1
fi

########################################
# 6. Build BLAST db
########################################
echo "[INFO] Building BLAST database..."
makeblastdb -in "$CHRM_FA" -dbtype nucl >/dev/null

########################################
# 7. BLAST merged loci to chrM
########################################
echo "[INFO] Running blastn..."
blastn \
  -query "$MERGED_FA" \
  -db "$CHRM_FA" \
  -out "$BLAST_OUT" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

########################################
# 8. Filter BLAST hits
########################################
echo "[INFO] Filtering BLAST hits (length >= ${MIN_ALIGN_LEN}, pident >= ${MIN_PIDENT})..."
awk -v min_len="$MIN_ALIGN_LEN" -v min_pid="$MIN_PIDENT" \
  '$4 >= min_len && $3 >= min_pid' \
  "$BLAST_OUT" > "$FILTERED_OUT"

########################################
# 9. Select best hit per locus by bitscore
########################################
echo "[INFO] Selecting best hit per locus..."
if [[ -s "$FILTERED_OUT" ]]; then
  sort -k1,1 -k12,12nr "$FILTERED_OUT" | awk '!seen[$1]++' > "$BESTHIT_RAW"
else
  : > "$BESTHIT_RAW"
fi

########################################
# 10. Join merged summary with chrM best hit
########################################
echo "[INFO] Joining merged summary with chrM best hit..."
awk '
BEGIN{
  OFS="\t";
}
FNR==NR{
  if(FNR==1) next;
  locus=$9;
  summary[locus]=$0;
  next;
}
{
  best[$1]=$0;
}
END{
  for(locus in summary){
    split(summary[locus], a, "\t");
    if(locus in best){
      split(best[locus], b, "\t");

      s1=b[9]+0;
      s2=b[10]+0;
      hit_start=(s1<s2 ? s1 : s2);
      hit_end=(s1>s2 ? s1 : s2);
      strand=(s1<=s2 ? "+" : "-");

      print a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9], \
            b[2],b[3],b[4],b[5],b[6],b[7],b[8],b[9],b[10],hit_start,hit_end,strand,b[11],b[12];
    } else {
      print a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9], \
            ".",".",".",".",".",".",".",".",".",".",".",".",".",".";
    }
  }
}
' "$MERGED_SUMMARY" "$BESTHIT_RAW" | sort -k1,1 -k2,2n > "$TMP_BODY"

{
  echo -e "nuclear_chrom\tnuclear_start\tnuclear_end\tn_blocks\tsum_nreads\tmax_nreads\tmean_nreads\tmean_MAPQ\tlocus_id\tchrM_hit\tpident\talign_length\tmismatch\tgapopen\tqstart\tqend\tchrM_start\tchrM_end\tchrM_hit_start\tchrM_hit_end\tstrand\tevalue\tbitscore"
  cat "$TMP_BODY"
} > "$BESTHIT_FINAL"

rm -f "$TMP_BODY"

########################################
# Done
########################################
echo "[INFO] Done."
echo "[INFO] Final output:"
echo "  $BESTHIT_FINAL"
echo
echo "[INFO] Intermediate files:"
echo "  $SORTED_BED"
echo "  $MERGED_BED"
echo "  $MERGED_SUMMARY"
echo "  $MERGED_FA"
echo "  $CHRM_FA"
echo "  $BLAST_OUT"
echo "  $FILTERED_OUT"
echo "  $BESTHIT_RAW"