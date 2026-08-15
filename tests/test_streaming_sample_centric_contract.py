#!/usr/bin/env python3
"""Production invariants for the stable per-sample streaming launcher."""
from pathlib import Path
import re

SCRIPT = (Path(__file__).resolve().parents[1] / "launch_pipeline_streaming_per_sample.sh").read_text()

assert 'SAMPLE_WORK_ROOT="${NF_BASE_WORK_DIR}/${SAMPLE_ID}"' in SCRIPT
assert "wave_" not in SCRIPT
assert 'flock -n 9 || die "Another worker is already using ${SAMPLE_ID}"' in SCRIPT
assert '^[A-Za-z0-9][A-Za-z0-9._-]*$' in SCRIPT
assert 'case "$sample_real" in "$root_real"/*)' in SCRIPT

for name in ("sample_manifest_sha256", "pipeline_config_sha256", "pipeline_git_commit",
             "reference_fingerprint", "important_parameters_sha256", "combined_fingerprint"):
    assert name in SCRIPT
assert 'cmp -s "$fingerprint_tmp" "${METADATA_DIR}/fingerprint.tsv"' in SCRIPT
assert '${SAMPLE_ID}.stale.$(date -u +%Y%m%dT%H%M%SZ).$$' in SCRIPT
for line in ('INFO: sample_fingerprint_match=', 'INFO: resume_mode=', 'INFO: sample_work_root='):
    assert line in SCRIPT

assert "#SBATCH --requeue" in SCRIPT and "#SBATCH --signal=B:USR1@300" in SCRIPT
handler = SCRIPT[SCRIPT.index("handle_sample_requeue()") : SCRIPT.index("trap handle_sample_requeue USR1")]
for item in ("TIMEOUT_SIGNAL", "resume_eligible\\t1", 'kill -TERM "$ACTIVE_CHILD_PID"',
             "for _ in $(seq 1 60)", 'kill -KILL "$ACTIVE_CHILD_PID"',
             'wait "$ACTIVE_CHILD_PID"', 'write_running_marker REQUEUE_REQUESTED',
             'scontrol requeue "$element_id"'):
    assert item in handler
assert "safe_remove" not in handler

for stage in ("pre", "numt", "round1", "round2"):
    assert f'run_stage {stage} "${{LOG_DIR}}/{stage}.log"' in SCRIPT
    assert f'clean_stage "$' + stage.upper() + '_WORK_DIR"' in SCRIPT
assert '-resume -w "$2"' in SCRIPT

final_check = SCRIPT[SCRIPT.index("final_outputs_complete()") : SCRIPT.index("write_numt_config()")]
for item in ("CRAM_PATH", "CRAI_PATH", "samtools quickcheck", "ROUND2_VCF",
             'gzip -t "$ROUND2_VCF"', "ROUND2_COVERAGE", "numt_decoy_coverage_path", "ROUND2_MTCN"):
    assert item in final_check
success = SCRIPT[SCRIPT.index("if (( success )); then") : SCRIPT.index('log "Completed per-sample')]
assert success.index(".complete.tsv") < success.index('safe_remove_sample_work "$SAMPLE_ID"')

assert 'max_attempts=$((1 + IMMEDIATE_SAMPLE_RETRIES))' in SCRIPT
assert '(( REQUEUE_IN_PROGRESS == 0 )) || exit 0' in SCRIPT
assert 'log "Immediate retry ${attempt}/${IMMEDIATE_SAMPLE_RETRIES}"' in SCRIPT
for classification in ("MISSING_REFERENCE", "MALFORMED_METADATA", "UNSUPPORTED_REFERENCE",
                       "UNSAFE_PATH", "FINGERPRINT_GENERATION_FAILED"):
    assert classification in SCRIPT
for epoch in ("first_failure_epoch", "last_failure_epoch", "last_attempt_epoch"):
    assert epoch in SCRIPT
assert re.search(r'FIRST_FAILURE_EPOCH=.*awk.*first_failure_epoch', SCRIPT)

for function in ("safe_remove_sample_stage_work()", "safe_remove_sample_work()", "normalize_manifest()"):
    assert function in SCRIPT
assert 'rm -rf -- "$target"' in SCRIPT
assert 'sbatch_args=(--export=ALL "--array=1-${NUM_SAMPLES}%${MAX_CONCURRENT}")' in SCRIPT
assert 'submission=$(sbatch "${sbatch_args[@]}" "$0")' in SCRIPT
assert "printf '%s\\n' \"$submission\"" in SCRIPT
assert 'SAMPLE_LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$NORMALIZED_SAMPLE_LIST")' in SCRIPT

slurm_id = SCRIPT[SCRIPT.index("current_slurm_element_id()") : SCRIPT.index("write_running_marker()")]
assert "printf '%s_%s\\n'" in slurm_id
assert "printf '%s\\n' \"${SLURM_JOB_ID:-}\"" in slurm_id
assert 'worker_state\\tTERMINAL_FAILED' in SCRIPT
assert 'write_running_marker IMMEDIATE_RETRY "$attempt"' in SCRIPT
assert '[[ -z "$FAILURE_CLASS" || "$FAILURE_CLASS" == UNKNOWN ]]' in SCRIPT
assert 'mkdir -p "$(dirname "$stage_log")"' in SCRIPT
assert 'STREAM_SMOKE_TEST=1 requires exactly one normalized sample' in SCRIPT
print("sample-centric streaming source contract: PASS")

# Round1 workflow and wrapper delegate to one checker; the obsolete WDL VCF is
# intentionally not a stage-completion requirement.
ROUND1 = (Path(__file__).resolve().parents[1] / "primate_pipeline_numt_decoy_round1.nf").read_text()
CHECKER = (Path(__file__).resolve().parents[1] / "scripts/round1_outputs_complete.sh").read_text()
assert 'scripts/round1_outputs_complete.sh' in ROUND1
round1_function = SCRIPT[SCRIPT.index("round1_complete()") : SCRIPT.index("round2_complete()")]
assert 'scripts/round1_outputs_complete.sh' in round1_function
assert "round_1_variant_calling_decoy" not in round1_function
assert "numt_decoy.clean.final.split.vcf" not in CHECKER
assert 'FAILURE_CLASS=OUTPUT_INCOMPLETE' in SCRIPT
assert 'FAILURE_REASON="Authoritative Round 1 outputs incomplete"; FAILURE_NONRETRYABLE=1' in SCRIPT

# A non-zero Nextflow stage remains retryable, while deterministic successful
# execution/output-contract mismatches do not consume an immediate retry.
run_stage_body = SCRIPT[SCRIPT.index("run_stage()") : SCRIPT.index("nf()")]
assert "FAILURE_NONRETRYABLE=1" not in run_stage_body.split('if [[ "$stage" == pre')[0]
assert '(( FAILURE_NONRETRYABLE == 0 )) || break' in SCRIPT

# Other stage paths match their current publishers. Round2 workflow skip logic
# now checks the same three concrete outputs rather than directory existence.
ROUND2 = (Path(__file__).resolve().parents[1] / "primate_pipeline_round2_consensus_NUMT.nf").read_text()
for suffix in ("round2.original_coords.clean.final.split.vcf.gz",
               "round2.original_coords.per_base_coverage.tsv", "round2.mtcn.tsv"):
    assert suffix in SCRIPT and suffix in ROUND2
assert 'requiredRound2Paths.every { it.exists() && it.size() > 0 }' in ROUND2
