# Streaming resume input audit

The launcher keeps `NF_BASE_WORK_DIR/<sample>` stable and treats its outer
fingerprint as provenance only. Nextflow remains responsible for cache
validation at process granularity.

## Input coverage

* `preprocessing.nf` reads the sample manifest through a `path` channel and
  stages downloaded FASTQs as `path(reads)`. Chunk routing values and alignment
  parameters are interpolated into the relevant process scripts. The launcher
  reference-bundle fingerprint is also interpolated into every process that
  reads the whole-genome reference (`ALIGN_AND_SORT`, both chunked alignment
  processes, and `BAM_TO_CRAM`). This covers reference content/index changes at
  a stable path, plus `ENABLE_CHUNKED_ALIGNMENT`, `FASTQ_SIZE_THRESHOLD_GB`, and
  `READS_PER_CHUNK` (mapped by the production Nextflow configuration to the
  workflow parameters).
* `numt_detection/numt_end2end.nf` receives its generated NUMT configuration as
  a `path` input. That configuration contains the sample/reference mapping and
  output/reference locations. Its sole task additionally interpolates the
  launcher reference-bundle fingerprint, covering whole-genome, chrM, and
  nuclear-only reference files and indexes that can change in place.
* Round 1 stages receive CRAM/CRAI, NUMT BED, generated references, VCFs, and
  consensus files through `path` inputs once located. `PREPARE_DECOY_REFERENCE`,
  the stage that directly opens the external whole-genome reference, also
  interpolates the reference-bundle fingerprint. Its derived reference bundle
  is passed downstream as `path` inputs.
* Round 2 first materializes the Round 1 BAM, VCF/index, original/consensus NUMT
  FASTA/index, and NUMT VCF/index as declared `path` outputs, and downstream
  processes consume those files as `path` inputs. `BUILD_CONSENSUS_REFERENCE`
  directly opens chrM and therefore also interpolates the reference-bundle
  fingerprint. Calling/alignment parameters are interpolated into the task
  scripts that use them.

The workflow source itself is part of Nextflow's task hash. Consequently a
Round2-only source edit can invalidate Round 2 tasks without hiding the
independent preprocessing, NUMT, or Round 1 cache directories. Existing
published-output validation remains a separate launcher decision and was not
changed by the resume fix.
