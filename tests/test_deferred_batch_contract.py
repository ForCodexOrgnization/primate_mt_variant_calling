#!/usr/bin/env python3
"""Contracts for throughput-first terminal failure handling."""
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
chain = (ROOT / "launch_pipeline_all_per_batch.sh").read_text()
nf = (ROOT / "preprocessing.nf").read_text()

for function in ("write_batch_status_marker()", "write_batch_complete_marker()",
                 "write_batch_failure_marker()", "classify_stage_failure()"):
    assert function in chain
assert "mktemp" in chain and 'mv -f "${tmp}" "${target}"' in chain
assert "REQUEUE_REQUESTED" in chain and "TIMEOUT_SIGNAL" in chain
assert "OUTPUT_INCOMPLETE" in chain and "CRAM_VALIDATION_FAILED" in chain
assert "failed_samples.tsv" in chain
for classification in ("TRANSIENT_IO", "BAD_INPUT", "DETERMINISTIC", "UNKNOWN"):
    assert classification in chain
assert re.search(r'wait "\$\{ACTIVE_CHILD_PID\}"\s+local rc=\$\?', chain)
assert '> >(tee -a "${stage_log}")' in chain
assert "successful batch cache cleanup completed" in chain
assert "refusing unsafe batch cache cleanup" in chain

split = nf[nf.index("process SPLIT_FASTQ_CHUNKS"):nf.index("process ALIGN_AND_SORT_CHUNK")]
for diagnostic in ("split_failure_diagnostics.txt", "partial_chunks_sizes.tsv",
                   "df -h .", "df -ih .", "du -sh .", "quota -s", "ulimit -a"):
    assert diagnostic in split
assert 'trap - INT TERM HUP QUIT ERR' in split
assert 'rm -rf chunks chunks.tsv chunk_ids.txt' in split

# A regression test for the exact set -e/process-substitution shape used by
# run_child_stage: wait must expose the child status, never tee's status.
probe = r'''set -euo pipefail
run() {
  set +e
  bash -c 'echo child-output; exit 122' > >(tee /dev/null) 2> >(tee /dev/null >&2) &
  pid=$!
  wait "$pid"
  rc=$?
  set -e
  return "$rc"
}
if run; then exit 99; else test "$?" -eq 122; fi
'''
subprocess.run(["bash", "-c", probe], check=True, capture_output=True, text=True)
print("deferred batch source contract: PASS")
