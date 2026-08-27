#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
script="${repo}/launch_pipeline_streaming_per_sample.sh"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/bin" "$tmp/refs/global" "$tmp/refs/mt" "$tmp/refs/nuclear" "$tmp/work"
printf 'x\n' >"$tmp/config"
printf 'sample_id\treference_name\ns1\tref\n' >"$tmp/samples.tsv"

make_tool() {
    cat >"$tmp/bin/$1" <<EOF
#!/usr/bin/env bash
[[ '$1' == nextflow ]] && echo
echo '$1 test-version'
EOF
    chmod +x "$tmp/bin/$1"
}
make_tool nextflow
make_tool samtools
cat >"$tmp/bin/sbatch" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$*" >"$SBATCH_CAPTURE"
printf '%s\n' "$PIPELINE_REPO_DIR" >"${SBATCH_CAPTURE}.repo"
echo 'Submitted batch job 12345'
EOF
chmod +x "$tmp/bin/sbatch"

run_coordinator() {
    env PATH="$tmp/bin:/usr/bin:/bin" SBATCH_CAPTURE="$tmp/sbatch.args" \
        FULL_SAMPLE_LIST="$tmp/samples.tsv" PRE_OUTPUT_DIR="$tmp/pre" ROUND_OUTPUT_DIR="$tmp/round" \
        NF_BASE_WORK_DIR="$tmp/work" NF_CONFIG_FILE="$tmp/config" GLOBAL_REF_DIR="$tmp/refs/global" \
        REF_DIR="$tmp/refs/mt" NUCLEAR_ONLY_REF_DIR="$tmp/refs/nuclear" "$@" "$script"
}

output=$(run_coordinator 2>&1)
grep -Fx 'Submitted batch job 12345' <<<"$output" >/dev/null
grep -F -- '--array=1-1%10' "$tmp/sbatch.args" >/dev/null
grep -F 'INFO: nextflow_executable=' <<<"$output" >/dev/null
grep -F 'INFO: samtools_executable=' <<<"$output" >/dev/null
grep -F 'INFO: nextflow_version=nextflow test-version' <<<"$output" >/dev/null
grep -Fx "$repo" "$tmp/sbatch.args.repo" >/dev/null
grep -F "INFO: pipeline_repo_dir=$repo" <<<"$output" >/dev/null

# The run-manager can request an initially held array without replacing the
# existing per-array concurrency throttle.
output=$(run_coordinator STREAMING_SUBMIT_HELD=1 MAX_CONCURRENT=7 2>&1)
grep -F -- '--array=1-1%7' "$tmp/sbatch.args" >/dev/null
grep -F -- '--hold' "$tmp/sbatch.args" >/dev/null

output=$(run_coordinator STREAMING_SUBMIT_HELD=0 2>&1)
grep -F -- '--array=1-1%10' "$tmp/sbatch.args" >/dev/null
if grep -F -- '--hold' "$tmp/sbatch.args" >/dev/null; then
    echo 'normal streaming submission unexpectedly held' >&2; exit 1
fi

output=$(run_coordinator 2>&1)
grep -F -- '--array=1-1%10' "$tmp/sbatch.args" >/dev/null
if grep -F -- '--hold' "$tmp/sbatch.args" >/dev/null; then
    echo 'default streaming submission unexpectedly held' >&2; exit 1
fi

rm -f "$tmp/sbatch.args" "$tmp/sbatch.args.repo"
if run_coordinator STREAMING_SUBMIT_HELD=abc >"$tmp/invalid-held.out" 2>&1; then
    echo 'invalid STREAMING_SUBMIT_HELD unexpectedly succeeded' >&2; exit 1
fi
grep -F 'ERROR: STREAMING_SUBMIT_HELD must be 0 or 1' "$tmp/invalid-held.out" >/dev/null
[[ ! -e "$tmp/sbatch.args" ]] || { echo 'invalid STREAMING_SUBMIT_HELD invoked sbatch' >&2; exit 1; }

# A copied Slurm script must use the explicitly exported repository, while a
# spool copy without that stable value must fail before any submission.
spool_dir="/var/spool/slurmd/job-test-$$"
mkdir -p "$spool_dir"
cp "$script" "$spool_dir/launcher.sh"
output=$(run_coordinator PIPELINE_REPO_DIR="$repo" "$spool_dir/launcher.sh" 2>&1)
grep -F "INFO: pipeline_repo_dir=$repo" <<<"$output" >/dev/null
if run_coordinator PIPELINE_REPO_DIR= "$spool_dir/launcher.sh" >"$tmp/spool.out" 2>&1; then
    echo 'unstable spool repository unexpectedly succeeded' >&2; exit 1
fi
grep -F 'Refusing to use Slurm spool directory as pipeline repository' "$tmp/spool.out" >/dev/null
rm -rf "$spool_dir"

# Every workflow is validated before tools or sbatch are invoked.
fake_repo="$tmp/incomplete-repo"
mkdir -p "$fake_repo/numt_detection" "$fake_repo/scripts"
workflows=(preprocessing.nf numt_detection/numt_end2end.nf primate_pipeline_numt_decoy_round1.nf primate_pipeline_round2_consensus_NUMT.nf scripts/round1_outputs_complete.sh)
for workflow in "${workflows[@]}"; do
    printf 'workflow test\n' >"$fake_repo/$workflow"
done
for workflow in "${workflows[@]}"; do
    mv "$fake_repo/$workflow" "$fake_repo/$workflow.saved"
    if run_coordinator PIPELINE_REPO_DIR="$fake_repo" >"$tmp/missing-workflow.out" 2>&1; then
        echo "missing workflow unexpectedly succeeded: $workflow" >&2; exit 1
    fi
    grep -F "Required pipeline script is missing: $fake_repo/$workflow" "$tmp/missing-workflow.out" >/dev/null
    mv "$fake_repo/$workflow.saved" "$fake_repo/$workflow"
done

rm "$tmp/bin/nextflow"
if run_coordinator >"$tmp/missing.out" 2>&1; then
    echo 'missing nextflow unexpectedly succeeded' >&2; exit 1
fi
grep -F 'NEXTFLOW_MODULE is empty' "$tmp/missing.out" >/dev/null

module_tools="$tmp/module-tools"; mkdir -p "$module_tools"
make_tool nextflow; mv "$tmp/bin/nextflow" "$module_tools/nextflow"
export MODULE_TOOLS="$module_tools"
module() { [[ "$1" == load && "$2" == test-nextflow ]] && PATH="$MODULE_TOOLS:$PATH"; }
export -f module
output=$(run_coordinator NEXTFLOW_MODULE=test-nextflow 2>&1)
grep -F 'INFO: nextflow_executable=' <<<"$output" >/dev/null

# A configured module is loaded even when an older executable is already in
# PATH, and the executable/version are resolved again after the load.
module_samtools="$tmp/module-samtools"; mkdir -p "$module_samtools"
cat >"$module_samtools/samtools" <<'EOF'
#!/usr/bin/env bash
echo 'samtools 1.21'
EOF
chmod +x "$module_samtools/samtools"
export MODULE_SAMTOOLS="$module_samtools" MODULE_CAPTURE="$tmp/module.loads"
module() {
    [[ "$1" == load ]] || return 1
    printf '%s\n' "$2" >>"$MODULE_CAPTURE"
    case "$2" in
        test-nextflow) PATH="$MODULE_TOOLS:$PATH" ;;
        test-samtools) PATH="$MODULE_SAMTOOLS:$PATH" ;;
        *) return 1 ;;
    esac
}
export -f module
output=$(run_coordinator NEXTFLOW_MODULE=test-nextflow SAMTOOLS_MODULE=test-samtools 2>&1)
grep -Fx 'test-samtools' "$MODULE_CAPTURE" >/dev/null
grep -F "INFO: samtools_executable=$module_samtools/samtools" <<<"$output" >/dev/null
grep -F 'INFO: samtools_version=samtools 1.21' <<<"$output" >/dev/null
grep -F 'hash -r' "$script" >/dev/null

printf 's2\tref\n' >>"$tmp/samples.tsv"
if run_coordinator NEXTFLOW_MODULE=test-nextflow STREAM_SMOKE_TEST=1 >"$tmp/smoke.out" 2>&1; then
    echo 'multi-sample smoke test unexpectedly succeeded' >&2; exit 1
fi
grep -F 'requires exactly one normalized sample' "$tmp/smoke.out" >/dev/null

echo 'streaming launcher integration tests: PASS'
