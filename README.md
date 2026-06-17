# primate_mt_variant_calling

## mtCN output

The consensus NUMT round 2 pipeline calculates mitochondrial copy number (mtCN) after round 2 chrM-assigned BAMs are generated. The result is written to:

```text
<outdir>/<sample_id>/round_2/mtcn/<sample_id>.round2.mtcn.tsv
```

The mtCN file contains round 2 mitochondrial mean coverage, WGS nuclear mean coverage, and `mtcn`, using:

```text
mtCN = 2 * mt_mean_coverage / nuclear_mean_coverage
```

Mitochondrial coverage is computed from the round 2 standard chrM-assigned BAM. Nuclear coverage is computed with `mosdepth` on the original WGS CRAM over all non-mitochondrial contigs in the whole-genome reference FASTA index. The mitochondrial contig is selected with `--mt_contig` (default: `chrM`).
