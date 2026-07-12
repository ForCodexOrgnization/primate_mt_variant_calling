#!/usr/bin/env bash
set -euo pipefail

# Deprecated. NUMT environment is now managed by nextflow.config / nextflow_mcc.config.
# This file is intentionally a no-op and is kept only for backward-compatible
# manual invocations that may still source it. Configure module loads and conda
# activation through params.numt_env_before_script and the numt_related label.

return 0 2>/dev/null || exit 0
