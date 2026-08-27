# Per-run FASTQ download deployment

The former preprocessing graph resolved and downloaded every run in one
sample-level `DOWNLOAD_FASTQ` task. A retry therefore received a new work hash
and repeated early downloads. The new graph first completes one authoritative
`RESOLVE_READ_RUNS` manifest, fans it out to one independently cacheable
`DOWNLOAD_FASTQ_RUN` task per run, preserves standard/chunked run alignment,
and gates the sample merge on the complete expected run identity set.

## Safe rollout

1. Merge and validate this change without deploying it into a checkout used by
   active sample workers.
2. Quiesce submissions, or create/update a versioned stable checkout at the
   merged commit. Never change pipeline commits between stages of one sample.
3. Direct only newly submitted or deliberately restarted samples to that stable
   checkout. Let active old-code samples finish, or cancel and restart them as a
   conscious operator action using the new commit.
4. Confirm published CRAM, CRAI, and completion markers with the existing
   validator before allowing downstream stages or cleanup.
5. Old failed monolithic `DOWNLOAD_FASTQ` work directories (including partial
   SRS003155 attempts) are not compatible cache entries. They may be cleaned
   manually only after the new deployment and outputs are validated. This
   pipeline change never deletes them automatically.

The principal operational change is increased task count. Network concurrency
is bounded by `params.download_run_max_forks` (default `2` in both production
configs); scheduler capacity and ENA behavior should be monitored during rollout.
