# Chunk-level interruption/resume acceptance test

This is an integration test for a Slurm test partition with Nextflow, BWA, samtools,
and pigz modules available. It intentionally retains the work directory.

1. Make one fixed PE sample TSV and set `fastq_chunk_reads` low enough to produce at
   least three chunks. Use a dedicated, persistent `WORK=/scratch/.../resume-test`
   and set `withName:ALIGN_AND_SORT_CHUNK.maxForks = 1` in a test config.
2. Start the exact command below and wait until trace/timeline shows two successful
   `ALIGN_AND_SORT_CHUNK` tasks. Terminate **only the Nextflow coordinator** (or send
   its Slurm coordinator job `USR1`):

   ```bash
   nextflow run preprocessing.nf -c nextflow.config -profile cluster \
     -w "$WORK" -resume --sample_tsv "$BATCH" --outdir "$OUT" \
     --enable_chunked_alignment true --chunked_alignment_size_threshold_gb 0 \
     --fastq_chunk_reads 1000
   ```

3. Repeat the byte-for-byte identical command with the same `WORK` and `BATCH`.
   The log must show the split and first two alignments as `CACHED`; only unfinished
   chunks may be submitted. Save `-with-trace trace.resume.txt` output as evidence.
4. Inject a one-time nonzero exit into one chunk task in a test-only config/script;
   verify only that task retries. Likewise inject a merge failure and verify every
   chunk is cached on resume. Do not cancel Slurm child jobs manually: after the
   coordinator restart, Nextflow must reconcile them or safely rerun only tasks
   without a completed cache entry.
5. Compare new and deprecated streaming outputs with `samtools view -c`,
   `samtools flagstat`, and read-name/pair counts. Verify the run BAM and final CRAM
   with `samtools quickcheck`; verify `chunks.tsv` row count equals every row's
   `expected_chunks`.

The dependency-free `python3 tests/test_chunk_resume_contract.py` test guards the
wiring, stable naming, validation, task-local cleanup, and launcher resume contract.
The full interruption test cannot be faithfully emulated without a Slurm executor.
