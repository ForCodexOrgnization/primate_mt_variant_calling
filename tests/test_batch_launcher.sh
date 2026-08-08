#!/usr/bin/env bash
set -euo pipefail

repo=$(realpath -m "$(dirname "${BASH_SOURCE[0]}")/..")
script="${repo}/launch_pipeline_all_per_batch.sh"
tmp=$(mktemp -d)
spool_dir="/var/spool/slurmd/job-batch-test-$$"
trap 'rm -rf "$tmp" "$spool_dir"' EXIT
mkdir -p "$tmp/bin" "$tmp/work" "$tmp/log"

cat >"$tmp/bin/sbatch" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$*" >"$SBATCH_CAPTURE"
printf '%s\n' "$PIPELINE_REPO_DIR" >"${SBATCH_CAPTURE}.repo"
echo 'Submitted batch job 12345'
EOF
chmod +x "$tmp/bin/sbatch"

for i in $(seq -w 1 100); do printf 'sample%s\tref\n' "$i"; done >"$tmp/samples.tsv"

run_parent() {
    env PATH="$tmp/bin:/usr/bin:/bin" SBATCH_CAPTURE="$tmp/sbatch.args" \
        FULL_SAMPLE_LIST="$tmp/samples.tsv" NF_BASE_WORK_DIR="$tmp/work" \
        BATCH_LIST_DIR="$tmp/work/batches" LOG_DIR="$tmp/log" BATCH_SIZE=5 \
        CHAIN_CONCURRENT_BATCHES=3 "$@"
}

# Direct invocation resolves the checked-out repository, exports it to sbatch,
# creates every fixed-size batch, and submits exactly one bounded array.
output=$(run_parent "$script" 2>&1)
grep -F "INFO: pipeline_repo_dir=$repo" <<<"$output" >/dev/null
[[ $(find "$tmp/work/batches" -maxdepth 1 -type f -name 'sample_batch_*' | wc -l) -eq 20 ]]
grep -F -- "--array=0-19%3 $script" "$tmp/sbatch.args" >/dev/null
grep -Fx "$repo" "$tmp/sbatch.args.repo" >/dev/null

# An explicit stable path wins for a script copied into the Slurm spool, and
# represents the value inherited by a nested array worker.
mkdir -p "$spool_dir"
cp "$script" "$spool_dir/launcher.sh"
output=$(run_parent PIPELINE_REPO_DIR="$repo" "$spool_dir/launcher.sh" 2>&1)
grep -F "INFO: pipeline_repo_dir=$repo" <<<"$output" >/dev/null
grep -Fx "$repo" "$tmp/sbatch.args.repo" >/dev/null
grep -F -- " $repo/launch_pipeline_all_per_batch.sh" "$tmp/sbatch.args" >/dev/null

# Simulate a nested worker receiving the exported value.  An out-of-range task
# stops before run_chain, so this cannot execute a scientific pipeline.
if run_parent PIPELINE_REPO_DIR="$repo" SLURM_ARRAY_TASK_ID=999 "$spool_dir/launcher.sh" >"$tmp/worker.out" 2>&1; then
    echo 'out-of-range array worker unexpectedly succeeded' >&2
    exit 1
fi
grep -F "INFO: pipeline_repo_dir=$repo" "$tmp/worker.out" >/dev/null
grep -F 'Could not find batch file for task ID 999' "$tmp/worker.out" >/dev/null

# Without the inherited stable path, a spool copy is rejected before sbatch or
# any scientific launcher can run.
rm -f "$tmp/sbatch.args"
if run_parent PIPELINE_REPO_DIR= "$spool_dir/launcher.sh" >"$tmp/spool.out" 2>&1; then
    echo 'unstable spool repository unexpectedly succeeded' >&2
    exit 1
fi
grep -F 'Refusing to use Slurm spool directory as pipeline repository' "$tmp/spool.out" >/dev/null
[[ ! -e "$tmp/sbatch.args" ]]

# Required launchers/workflows are resolved and validated under the explicit
# repository before submission.
fake_repo="$tmp/incomplete-repo"
mkdir -p "$fake_repo/numt_detection"
if run_parent PIPELINE_REPO_DIR="$fake_repo" "$script" >"$tmp/missing.out" 2>&1; then
    echo 'incomplete repository unexpectedly succeeded' >&2
    exit 1
fi
grep -F "Required pipeline script is missing: $fake_repo/launch_pipeline_all_per_batch.sh" "$tmp/missing.out" >/dev/null
[[ ! -e "$tmp/sbatch.args" ]]

echo 'batch launcher integration tests: PASS'
