#!/usr/bin/env bash
# Submit a 10-sample x (30,50,100)x(101,202) HOM experiment on McCleary.
# One array index = one manifest row. Dependencies are per matching array task.
# Example: MAX_CONCURRENT=4 bash scripts/submit_hom_downsampling_batch.sh
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CALL_REPO="$(cd "$SCRIPT_DIR/.." && pwd)"
ANALYSIS_ROOT="${HOM_ANALYSIS_ROOT:-/vast/palmer/scratch/lake_nicole/lt692/primate_mtDNA_analysis_coverage_only}"
BASE="$ANALYSIS_ROOT/results/qc_coverage_only/downsampling_pilot"
MANIFEST="$BASE/candidate_downsampling_jobs.tsv"
PARTITION="${HOM_SLURM_PARTITION:-ycga}"
MAX_CONCURRENT="${MAX_CONCURRENT:-4}"
NEXTFLOW_MODULE="${NEXTFLOW_MODULE:-Nextflow/24.04.2}"
[[ "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid MAX_CONCURRENT" >&2; exit 2; }
[[ -s "$MANIFEST" ]] || { echo "ERROR: Missing manifest: $MANIFEST" >&2; exit 1; }
[[ -s "$ANALYSIS_ROOT/qc_analysis/scripts/run_candidate_downsampling_array.slurm" ]] || { echo "ERROR: Analysis downsampling script missing" >&2; exit 1; }
[[ -s "$CALL_REPO/scripts/run_hom_downsampling_batch_stage.sh" ]] || { echo "ERROR: Round1/Round2 array runner missing" >&2; exit 1; }
read -r SAMPLE_COUNT N < <(
python3 - "$MANIFEST" <<'PY'
import csv,sys
from collections import Counter
with open(sys.argv[1],newline='') as f:
    r=list(csv.DictReader(f, delimiter='\t'))
expected={(d,s) for d in (30,50,100) for s in (101,202)}
groups={}
for x in r:
    sample=x['sample']
    key=(int(x['target_depth']),int(x['seed']))
    if key in groups.setdefault(sample,set()):
        raise SystemExit('ERROR: Duplicate task: '+str((sample,key)))
    groups[sample].add(key)
if len(groups)!=10 or len(r)!=60 or any(v!=expected for v in groups.values()):
    raise SystemExit(f'ERROR: Manifest is not 10 samples x 3 depths x 2 seeds: {len(groups)} samples, {len(r)} tasks')
print(len(groups),len(r))
PY
)
mkdir -p "$BASE/logs"
echo "Batch: $SAMPLE_COUNT samples; $N tasks; partition=$PARTITION max_parallel_per_stage=$MAX_CONCURRENT"
echo "Downsampling -> Round1 -> Round2, task-index-matched aftercorr dependencies."
if [[ "${1:-}" == "--dry-run" ]]; then
    echo "Dry run only: no Slurm jobs submitted."
    exit 0
fi
[[ $# -eq 0 ]] || { echo "Usage: MAX_CONCURRENT=4 bash $0 [--dry-run]" >&2; exit 2; }
cd "$ANALYSIS_ROOT"
D_JOB="$(sbatch --parsable --partition="$PARTITION" --array="1-${N}%${MAX_CONCURRENT}" \
    --output="$BASE/logs/hom_ds_%A_%a.out" --error="$BASE/logs/hom_ds_%A_%a.err" \
    --export=ALL "$ANALYSIS_ROOT/qc_analysis/scripts/run_candidate_downsampling_array.slurm")"
D_JOB="${D_JOB%%;*}"
echo "Submitted downsampling array: $D_JOB"
R1_JOB="$(sbatch --parsable --partition="$PARTITION" --array="1-${N}%${MAX_CONCURRENT}" \
    --dependency="aftercorr:$D_JOB" --job-name=hom_batch_r1 \
    --cpus-per-task=2 --mem=8G --time=24:00:00 \
    --output="$BASE/logs/hom_batch_r1_%A_%a.out" --error="$BASE/logs/hom_batch_r1_%A_%a.err" \
    --export=ALL,HOM_CALL_REPO="$CALL_REPO",HOM_ANALYSIS_ROOT="$ANALYSIS_ROOT",NEXTFLOW_MODULE="$NEXTFLOW_MODULE" \
    "$CALL_REPO/scripts/run_hom_downsampling_batch_stage.sh" round1)"
R1_JOB="${R1_JOB%%;*}"
echo "Submitted Round1 array: $R1_JOB (aftercorr:$D_JOB)"
R2_JOB="$(sbatch --parsable --partition="$PARTITION" --array="1-${N}%${MAX_CONCURRENT}" \
    --dependency="aftercorr:$R1_JOB" --job-name=hom_batch_r2 \
    --cpus-per-task=2 --mem=8G --time=24:00:00 \
    --output="$BASE/logs/hom_batch_r2_%A_%a.out" --error="$BASE/logs/hom_batch_r2_%A_%a.err" \
    --export=ALL,HOM_CALL_REPO="$CALL_REPO",HOM_ANALYSIS_ROOT="$ANALYSIS_ROOT",NEXTFLOW_MODULE="$NEXTFLOW_MODULE" \
    "$CALL_REPO/scripts/run_hom_downsampling_batch_stage.sh" round2)"
R2_JOB="${R2_JOB%%;*}"
echo "Submitted Round2 array: $R2_JOB (aftercorr:$R1_JOB)"
printf 'stage\tjob_id\n' > "$BASE/hom_batch_job_ids.tsv"
printf 'downsampling\t%s\nround1\t%s\nround2\t%s\n' "$D_JOB" "$R1_JOB" "$R2_JOB" >> "$BASE/hom_batch_job_ids.tsv"
echo "Job IDs: $BASE/hom_batch_job_ids.tsv"
echo "Monitor: squeue -u \$USER -o '%.18i %.12P %.25j %.2t %.12M %R'"
