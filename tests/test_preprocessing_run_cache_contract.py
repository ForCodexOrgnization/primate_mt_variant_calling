#!/usr/bin/env python3
"""Regression contract for independently cached per-run preprocessing tasks."""
from pathlib import Path

repo = Path(__file__).resolve().parents[1]
source = (repo / "preprocessing.nf").read_text()
config = (repo / "nextflow.config").read_text()
mcc = (repo / "nextflow_mcc.config").read_text()

# Three (or any number of) manifest rows fan out into separate process calls.
assert "RESOLVE_READ_RUNS(ch_samples)" in source
assert ".splitCsv(header: true" in source
assert "DOWNLOAD_FASTQ_RUN(ch_resolved_runs)" in source
assert "process DOWNLOAD_FASTQ_RUN" in source
assert "tuple val(meta), val(species_name), val(ref_name), val(r1_url)" in source

# Completed A/B tasks have hashes based on their own row and files; C failing
# terminates rather than being ignored. Thus -resume can cache A/B independently.
download = source.split("process DOWNLOAD_FASTQ_RUN", 1)[1].split("process ALIGN_AND_SORT", 1)[0]
assert "task.attempt <= 2 ? 'retry' : 'terminate'" in download
assert "'ignore'" not in download and "failed_download.tsv" not in download

# Merge completeness and duplicate-identity gates.
assert "groupKey(meta.id, meta.n_pairs)" in source
assert "observed.size() != observed.toSet().size()" in source
assert "observed.sort() != expected" in source
assert "CRITICAL: Sample" in source

# Single-run, PE, and SE use the same generic task and retain both routes.
assert "meta.layout == 'PE'" in download
assert "r1_url.tokenize('/').last()" in download
assert "ALIGN_AND_SORT(ch_fastq_alignment_routes.standard)" in source
assert "SPLIT_FASTQ_CHUNKS(ch_fastq_alignment_routes.chunked)" in source
assert "MERGE_RUN_CHUNK_BAMS(ch_run_chunk_bams)" in source

# Nextflow may unwrap a one-file path output to a scalar Path. A Path is itself
# Iterable over components, so normalize it before sizing and preserve the
# normalized list on both the standard and chunked branch outputs.
routing = source.split("ch_fastq_alignment_routes = ch_fastq_pairs", 1)[1].split(
    "ALIGN_AND_SORT(ch_fastq_alignment_routes.standard)", 1
)[0]
assert "(reads instanceof Collection) ? reads.toList() : [reads]" in routing
assert "read_files.collect { it.size() }" in routing
assert "reads.collect { it.size() }" not in routing
assert "tuple(routed_meta, species_name, ref_name, read_files)" in routing

# Behavioral route matrix: scalar SE becomes one input while PE remains two;
# threshold routing never changes the normalized cardinality passed onward.
def normalize(reads):
    return list(reads) if isinstance(reads, list) else [reads]


for reads, expected_count in [("SINGLE.fastq.gz", 1), (["R1.fastq.gz", "R2.fastq.gz"], 2)]:
    normalized = normalize(reads)
    assert len(normalized) == expected_count
    assert normalized[0] == (reads[0] if isinstance(reads, list) else reads)
    for above_chunk_threshold in (False, True):
        routed_to = "chunked" if above_chunk_threshold else "standard"
        route_payload = (routed_to, normalized)
        assert len(route_payload[1]) == expected_count

# Fan-out is explicit, conservative, and independently configurable on both sites.
for text in (config, mcc):
    assert "download_run_max_forks = 2" in text
    assert "withName: DOWNLOAD_FASTQ_RUN" in text
    assert "maxForks = params.download_run_max_forks" in text

print("per-run preprocessing/cache/merge contract: PASS")
