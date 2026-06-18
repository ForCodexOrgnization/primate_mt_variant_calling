# primate_mt_variant_calling

## mtCN output

The consensus NUMT round 2 pipeline calculates mitochondrial copy number (mtCN) after round 2 chrM-assigned BAMs are generated. The result is written to:

```text
<outdir>/<sample_id>/round_2/mtcn/<sample_id>.round2.mtcn.tsv
```

The mtCN file contains both mean-based and median-based round 2 mtCN values:

```text
mtcn_mean = 2 * mt_mean_coverage / nuclear_mean_coverage
mtcn_median = 2 * mt_median_coverage / nuclear_median_coverage
```

Mitochondrial coverage is computed from the round 2 standard chrM-assigned BAM. Nuclear coverage is computed with `mosdepth` on the original WGS CRAM over all non-mitochondrial contigs in the whole-genome reference FASTA index. Mean nuclear coverage is length-weighted across non-mt contigs; median nuclear coverage is the length-weighted median of the non-mt contig coverage values reported by `mosdepth --by`. The mitochondrial contig is selected with `--mt_contig` (default: `chrM`).
