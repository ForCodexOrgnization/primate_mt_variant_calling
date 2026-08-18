#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
bash "$repo/tests/test_fastq_classification.sh"
bash "$repo/tests/test_read_run_resolver.sh"
bash "$repo/tests/test_round1_completion.sh"
python3 "$repo/tests/test_streaming_sample_centric_contract.py"
python3 "$repo/tests/test_round_contract_source.py"
python3 "$repo/tests/test_mtcn_nuclear_metrics.py"
bash "$repo/tests/test_streaming_launcher.sh"
bash "$repo/tests/test_batch_launcher.sh"
