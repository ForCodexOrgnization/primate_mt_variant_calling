#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage:
  bash run_numt_end2end.sh --config numt_pipeline.config

Run NUMT discovery and downstream best-hit analysis in one command.
All parameters and paths are read from the config file.
USAGE
}

CONFIG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$CONFIG" ]] || { echo "ERROR: --config required" >&2; exit 1; }
[[ -s "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }

# shellcheck disable=SC1090
source "$CONFIG"

: "${SAMPLE:=${SAMPLE_ID:-}}"
: "${SAMPLE:?missing SAMPLE in config}"
: "${DISCOVERY_OUTDIR:?missing DISCOVERY_OUTDIR in config}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BED_PATH="${DISCOVERY_OUTDIR}/${SAMPLE}.numt_candidates.bed"
BESTHIT_DIR="${BESTHIT_OUTDIR:-/nfs/roberts/pi/pi_njl27/lt692/primate_results_numt_besthit}"
HIGHCONF_PATH="${BESTHIT_DIR}/${SAMPLE}.highconf_numt.bed"

if [[ -e "$HIGHCONF_PATH" ]]; then
  echo "[$(date)] highconf output already exists, skipping sample: $HIGHCONF_PATH"
  exit 0
fi

DISCOVERY_ARGS=(--config "$CONFIG")
if [[ -e "$BED_PATH" && ! -s "$BED_PATH" ]]; then
  echo "[$(date)] Empty candidate BED detected without highconf output; forcing discovery rerun from scratch."
  DISCOVERY_ARGS+=(--force-rerun)
fi

echo "[$(date)] Running discovery..."
bash "${SCRIPT_DIR}/run_numt_discovery.sh" "${DISCOVERY_ARGS[@]}"

echo "[$(date)] Running best-hit analysis..."
bash "${SCRIPT_DIR}/process_numt_candidates_one_ready.sh" \
  --bed "$BED_PATH" \
  --config "$CONFIG"

echo "[$(date)] End-to-end pipeline finished."
