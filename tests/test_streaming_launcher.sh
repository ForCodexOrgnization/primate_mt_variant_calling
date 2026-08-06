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
echo '$1 test-version'
EOF
    chmod +x "$tmp/bin/$1"
}
make_tool nextflow
make_tool samtools
cat >"$tmp/bin/sbatch" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$*" >"$SBATCH_CAPTURE"
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

printf 's2\tref\n' >>"$tmp/samples.tsv"
if run_coordinator NEXTFLOW_MODULE=test-nextflow STREAM_SMOKE_TEST=1 >"$tmp/smoke.out" 2>&1; then
    echo 'multi-sample smoke test unexpectedly succeeded' >&2; exit 1
fi
grep -F 'requires exactly one normalized sample' "$tmp/smoke.out" >/dev/null

echo 'streaming launcher integration tests: PASS'
