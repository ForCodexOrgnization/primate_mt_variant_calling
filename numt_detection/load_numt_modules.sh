#!/usr/bin/env bash
set -euo pipefail

# Centralized environment setup for NUMT pipeline on HPC.
# Override these in parent shell before sourcing if needed.
: "${NUMT_MODULE_SAMTOOLS:=SAMtools/1.21-GCC-13.3.0}"
: "${NUMT_MODULE_BWA:=BWA/0.7.18-GCCcore-13.3.0}"
: "${NUMT_MODULE_BEDTOOLS:=BEDTools/2.31.1-GCC-13.3.0}"
: "${NUMT_MODULE_MINICONDA:=miniconda/24.11.3}"
: "${NUMT_CONDA_ENV:=blast_env}"

if command -v module >/dev/null 2>&1; then
  module load "$NUMT_MODULE_SAMTOOLS"
  module load "$NUMT_MODULE_BWA"
  module load "$NUMT_MODULE_BEDTOOLS"
  module load "$NUMT_MODULE_MINICONDA"
fi

if command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  conda activate "$NUMT_CONDA_ENV" || true
fi
