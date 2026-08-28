# primate_mt_variant_calling

## Sequencing-platform input policy

The preprocessing pipeline accepts actual paired-end FASTQ runs from an
explicit short-read platform whitelist: **ILLUMINA**, **DNBSEQ**, and
**BGISEQ**. The downloadable FASTQ structure determined by
`scripts/classify_ena_fastq.sh` is authoritative; ENA `library_layout` is not
used to admit single-end data. Actual SE runs and PacBio, Oxford Nanopore,
LS454, or other platforms are excluded with their original ENA provenance.

Each biological sample uses exactly one platform. All eligible PE runs from
the highest available priority are selected in the order ILLUMINA > DNBSEQ >
BGISEQ; otherwise-eligible lower-priority runs are recorded as excluded and
are never merged into the sample CRAM. SAM read groups retain `PL:ILLUMINA`
for Illumina and use the canonical `PL:DNBSEQ` for both DNBSEQ and BGISEQ.
The used-run provenance continues to record the original ENA platform and
instrument model, so BGISEQ remains identifiable as BGISEQ there.

## mtCN output

The consensus NUMT round 2 pipeline calculates mitochondrial copy number (mtCN) after round 2 chrM-assigned BAMs are generated. The result is written to:

```text
<outdir>/<sample_id>/round_2/mtcn/<sample_id>.round2.mtcn.tsv
```

The mtCN file contains both mean-based and median-based round 2 mtCN values. For compatibility with downstream summary tables, it includes alias columns named `mean_mt_coverage`, `mean_nuclear_coverage`, and `mean_mtCN` in addition to the existing `mt_mean_coverage`, `nuclear_mean_coverage`, and `mtcn_mean` columns:

```text
mtcn_mean = 2 * mt_mean_coverage / nuclear_mean_coverage
mtcn_median = 2 * mt_median_coverage / nuclear_median_coverage
```

Mitochondrial coverage is computed from the round 2 standard chrM-assigned BAM. Nuclear coverage is computed with `mosdepth` on the original WGS CRAM over all non-mitochondrial contigs in the whole-genome reference FASTA index. Mean nuclear coverage is length-weighted across non-mt contigs; median nuclear coverage is the length-weighted median of the non-mt contig coverage values reported by `mosdepth --by`. The mitochondrial contig is selected with `--mt_contig` (default: `chrM`).

## Round 1 NUMT consensus filtering

NUMT consensus construction follows a mtSwirl-like filtering strategy. NUMT variants are first called with GATK HaplotypeCaller on the decoy NUMT interval list, then SNPs and INDELs are hard-filtered separately. Heterozygous genotypes, homozygous-reference genotypes, and genotypes with DP below the configured lower bound are excluded. Only PASS non-reference NUMT variants from the left-aligned, split VCF are used to build the consensus NUMT FASTA.

The default genotype depth lower bound is 10 and can be configured when launching Nextflow:

```bash
--hc_dp_lower_bound 10
```

## Round 2 chrM-assigned BAM ambiguity filtering

Round 2 keeps the existing mtSwirl-like preprocessing: FASTQ is regenerated from the round 1 candidate BAM, reads are competitively realigned to the standard and shifted chrM+NUMT self-references, duplicates are marked, and primary chrM alignments are written to chrM-only BAMs.

The QNAME-level ambiguity filter for those final chrM-assigned BAMs is controlled with:

```bash
--round2_chrm_assignment_filter 1
```

Available settings:

1. `1` (default): enable filtering. Remove a whole QNAME if any retained chrM primary record has `XS > AS`; remove `AS > XS` records and their QNAME mates when the mate is missing from the retained chrM-primary set or the mate lacks either `AS` or `XS`; and remove partial `XS == AS` ambiguous records and their mates when the mate is missing, lacks either `AS` or `XS`, or has `XS >= AS`.
2. `2`: disable this filter and pass retained primary chrM alignments through to the final chrM-assigned BAM.

Each branch writes `*.chrM_assignment_filter.metrics.tsv`, `*.chrM_assignment_filter.drop_qnames.txt`, and `*.chrM_assignment_filter.drop_reasons.tsv` next to the round 2 alignment outputs. When filtering is disabled, the drop files are still emitted but contain no dropped QNAMEs.
