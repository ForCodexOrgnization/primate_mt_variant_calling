#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
summary="$tmp/metadata/used_short_read_runs.tsv"
update="$repo/scripts/update_used_short_read_runs.sh"

"$update" "$summary" SAMPLE Species Ref 'SRR3,SRR1,SRR1,SRR2'
awk -F '\t' 'NR==2 { exit !($1=="SAMPLE" && $4==3 && $5=="SRR1,SRR2,SRR3") } END { exit !(NR==2) }' "$summary"
"$update" "$summary" SAMPLE NewSpecies NewRef 'SRR9'
awk -F '\t' 'NR==2 { exit !($1=="SAMPLE" && $2=="NewSpecies" && $4==1 && $5=="SRR9") } END { exit !(NR==2) }' "$summary"

for n in $(seq 1 20); do "$update" "$summary" "S$n" Species Ref "SRR$n" & done
wait
[[ $(tail -n +2 "$summary" | wc -l) -eq 21 ]]
[[ $(cut -f1 "$summary" | tail -n +2 | sort -u | wc -l) -eq 21 ]]
awk -F '\t' 'NR==1 { if (NF != 5) exit 1; next } NF != 5 { exit 1 }' "$summary"

# Provenance is emitted in BAM_TO_CRAM only after final CRAM validation, and
# derives from the exact expected set carried through MERGE_BAMS.
python3 - "$repo/preprocessing.nf" <<'PY'
import pathlib, sys
s = pathlib.Path(sys.argv[1]).read_text()
p = s[s.index('process BAM_TO_CRAM'):]
assert p.index('samtools idxstats') < p.index('used_short_read_runs.tsv.tmp')
assert 'def expectedRunIds = meta.expected_run_ids.sort().join' in p
assert 'scripts/update_used_short_read_runs.sh' in p
assert p.index('cram.complete.tmp" "${meta.id}.cram.complete') < p.index('scripts/update_used_short_read_runs.sh')
PY
echo 'used short-read provenance: PASS'
