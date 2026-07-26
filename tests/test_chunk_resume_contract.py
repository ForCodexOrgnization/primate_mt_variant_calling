#!/usr/bin/env python3
"""Fast, dependency-free regression checks for the chunk-resume workflow contract."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
nf = (ROOT / "preprocessing.nf").read_text()
launch = (ROOT / "launch_pipeline_pre.sh").read_text()

workflow = nf[nf.index("workflow {"):nf.index("/*\n================================================================================\n    PROCESSES")]
assert "ALIGN_LARGE_FASTQ_STREAMING_CHUNKS(" not in workflow
for call in ("SPLIT_FASTQ_CHUNKS(ch_fastq_alignment_routes.chunked)",
             "ALIGN_AND_SORT_CHUNK(ch_alignment_chunks)",
             "MERGE_RUN_CHUNK_BAMS(ch_run_chunk_bams)"):
    assert call in workflow, call
assert "groupKey([sample_id: meta.id, pair_id: meta.pair_id], meta.n_chunks)" in workflow
assert "bams.size() != meta.n_chunks" in workflow
assert "ALIGN_AND_SORT.out.bam.mix(MERGE_RUN_CHUNK_BAMS.out.bam)" in workflow

header = r"sample_id\\tspecies_name\\tref_name\\trun_id\\tlayout\\tchunk_id\\tr1\\tr2\\texpected_chunks\\texpected_pairs"
assert header in nf
assert "--numeric-suffixes=1" in nf and "-a 6" in nf
assert re.search(r'process ALIGN_AND_SORT_CHUNK.*?samtools index.*?samtools quickcheck', nf, re.S)
assert "tmp_sort_${meta.id}_${meta.pair_id}_${meta.chunk_id}" in nf
assert re.search(r'process MERGE_RUN_CHUNK_BAMS.*?bams\.size\(\).*?samtools merge.*?samtools index.*?samtools quickcheck.*?mv', nf, re.S)

for token in ("#SBATCH --requeue", "#SBATCH --signal=B:USR1@300", "trap handle_requeue USR1",
              'NEXTFLOW_PID=$!', '-resume -w "${WORK_DIR}"', "samtools quickcheck"):
    assert token in launch, token
print("chunk-resume source contract: PASS")
