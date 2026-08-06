#!/usr/bin/env python3
"""Regression checks for Slurm USR1 requeue handling in per-batch chains."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
chain = (ROOT / "launch_pipeline_all_per_batch.sh").read_text()
pre = (ROOT / "launch_pipeline_pre.sh").read_text()

assert "#SBATCH --time=24:00:00" in chain
assert "#SBATCH --requeue" in chain
assert "#SBATCH --signal=B:USR1@300" in chain

for token in ('ACTIVE_CHILD_PID=""', 'ACTIVE_STAGE=""', 'ACTIVE_BATCH_NAME=""',
              'ACTIVE_BATCH_FILE=""', 'ACTIVE_STAGE_LOG=""', 'REQUEUE_IN_PROGRESS=0'):
    assert token in chain, token

for fn in ("current_slurm_element_id()", "run_child_stage()", "handle_chain_requeue()"):
    assert fn in chain, fn
handler = chain[chain.index("handle_chain_requeue()") : chain.index("trap handle_chain_requeue USR1")]
for token in (
    'if [[ "${REQUEUE_IN_PROGRESS}" == 1 ]]',
    'log "INFO: ACTIVE_STAGE=${ACTIVE_STAGE:-none}"',
    'log "INFO: Preserving NF_BASE_WORK_DIR=${NF_BASE_WORK_DIR}"',
    'kill -TERM "${ACTIVE_CHILD_PID}"',
    'for _ in $(seq 1 60)',
    'sleep 2',
    'kill -KILL "${ACTIVE_CHILD_PID}"',
    'scontrol requeue "${element_id}"',
    'exit 0',
):
    assert token in handler, token
assert "cleanup_work_dir_if_requested" not in handler

for stage in ("pre", "numt", "round1", "round2"):
    assert re.search(rf"run_child_stage {stage} \"\$\{{LOG_DIR\}}/.*?\.{stage}\.log\"", chain), stage

pre_call = re.search(r"run_child_stage pre .*? env \\\n(?P<body>.*?)bash \"\$\{PRE_LAUNCH_SCRIPT\}\"", chain, re.S).group("body")
for token in (
    'BATCH_FILE="${batch_file}"',
    'BATCH_ID="${batch_name}"',
    'FULL_SAMPLE_LIST="${batch_file}"',
    'NF_BASE_WORK_DIR="${PRE_NF_BASE_WORK_DIR}"',
    'NF_CONFIG_FILE="${NF_CONFIG_FILE}"',
    'DEFER_WORK_DIR_CLEANUP=1',
    'CHAIN_MANAGED_REQUEUE=1',
):
    assert token in pre_call, token

numt = chain[chain.index("run_numt_nextflow()") : chain.index("validate_pre_to_round1()")]
assert "local numt_cmd=(" in numt
assert "run_child_stage numt \"${LOG_DIR}/${batch_name}.numt.log\" env" in numt
assert "rm -f" not in numt, "NUMT must not remove Nextflow cache metadata before retry"

for stage, script_var in (("round1", "ROUND1_LAUNCH_SCRIPT"), ("round2", "ROUND2_LAUNCH_SCRIPT")):
    start = chain.index(f"run_child_stage {stage} \"${{LOG_DIR}}/${{batch_name}}.{stage}.log\" env")
    end = chain.index(f'bash "${{{script_var}}}"', start)
    block = chain[start:end]
    assert "CHAIN_MANAGED_REQUEUE=1" in block

assert 'if [[ "${CHAIN_MANAGED_REQUEUE:-0}" != "1" ]]; then' in pre
assert "trap handle_requeue USR1" in pre
assert "INFO: Requeue is managed by launch_pipeline_all_per_batch.sh" in pre
batch_mode = pre[pre.index('if [ -n "${BATCH_FILE:-}" ]; then') : pre.index('if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then')]
for token in (
    'RUN_DIR="${BATCH_ID:-$(basename "${BATCH_FILE}")}"',
    'WORK_DIR="${NF_BASE_WORK_DIR}/${RUN_DIR}"',
    'INFO: CHAIN_MANAGED_REQUEUE=${CHAIN_MANAGED_REQUEUE:-0}',
    'INFO: Persistent Nextflow work directory=${WORK_DIR}',
    'INFO: Nextflow resume enabled',
    'INFO: BATCH_FILE=${ACTIVE_BATCH_FILE}',
    'INFO: BATCH_ID=${ACTIVE_BATCH_ID}',
):
    assert token in batch_mode, token
run_nextflow = pre[pre.index("run_nextflow()") : pre.index("verify_batch_outputs()")]
for token in ('local nf_cmd=(', '-resume', '-w "${WORK_DIR}"', '"${nf_cmd[@]}" &', 'NEXTFLOW_PID=$!'):
    assert token in run_nextflow, token
assert "rm -rf \"${WORK_DIR}\"" not in pre[pre.index('run_nextflow ||') : pre.index('echo "BATCH_FILE mode completed successfully')]
print("chain requeue source contract: PASS")
