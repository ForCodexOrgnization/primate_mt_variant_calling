#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mock="$tmp/curl"
cat >"$mock" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
out= url= accession=
args=("$@")
for ((i=0; i<${#args[@]}; i++)); do
  case ${args[i]} in
    -o) out=${args[++i]};;
    http*) url=${args[i]};;
    accession=*) accession=${args[i]#accession=};;
  esac
done
header=$'run_accession\tlibrary_layout\tfastq_ftp\tfastq_md5'
emit_run() {
  case $1 in
    SRR1746970) printf '%s\nSRR1746970\tPAIRED\th/SRR1746970_1.fastq.gz;h/SRR1746970_2.fastq.gz\tm1;m2\n' "$header";;
    SRR1746971) printf '%s\nSRR1746971\tPAIRED\th/SRR1746971.fastq.gz;h/SRR1746971_1.fastq.gz;h/SRR1746971_2.fastq.gz\tmo;m1;m2\n' "$header";;
    SRR1) printf '%s\nSRR1\tPAIRED\th/SRR1_1.fastq.gz;h/SRR1_2.fastq.gz\ta;b\n' "$header";;
    SRR999999) printf '%s\nSRR999999\tPAIRED\th/SRR999999_1.fastq.gz;h/SRR999999_2.fastq.gz\t\n' "$header";;
  esac
}
case $url in
  */filereport)
    if [[ ${MODE:-fallback} == direct_run && $accession == SRR1 ]]; then emit_run SRR1 >"$out"
    elif [[ ${MODE:-fallback} == direct_sample && $accession == ERS1 ]]; then emit_run SRR1 >"$out"
    elif [[ ${MODE:-fallback} == mixed_discovery && $accession == ERS_MIXED ]]; then
      { emit_run SRR1; printf '\tPAIRED\th/partial.fastq.gz\t\n'; } >"$out"
    elif [[ ${MODE:-fallback} =~ ^incomplete_ && $accession == ERS14600391 ]]; then printf '%s\n\tPAIRED\t\t\n' "$header" >"$out"
    elif [[ $accession =~ ^SRR ]]; then emit_run "$accession" >"$out"
    else printf '%s\n' "$header" >"$out"; fi;;
  */search)
    if [[ ${MODE:-fallback} =~ ^incomplete_ ]]; then printf '%s\nNOT_A_RUN\t\th/partial.fastq.gz\t\n' "$header" >"$out"
    else printf '%s\n' "$header" >"$out"; fi;;
  */esearch.fcgi)
    if [[ ${MODE:-fallback} == none ]]; then printf '{"esearchresult":{"idlist":[]}}' >"$out"
    else printf '{"esearchresult":{"idlist":["1"]}}' >"$out"; fi;;
  */efetch.fcgi)
    if [[ ${MODE:-fallback} == incomplete_strict ]]; then printf '<RUN_SET><RUN accession="SRR999999"/></RUN_SET>' >"$out"
    else printf '<RUN_SET><RUN accession="SRR1746971"/><RUN accession="SRR1746970"/><RUN accession="SRR1746971"/></RUN_SET>' >"$out"; fi;;
  *) exit 3;;
esac
MOCK
chmod +x "$mock"

CURL_BIN="$mock" "$repo/scripts/resolve_read_runs.sh" SAMN03275524 "$tmp/report"
[[ $(tail -n +2 "$tmp/report" | wc -l) -eq 2 ]]
[[ $(sed -n '2p' "$tmp/report" | cut -f1) == SRR1746970 ]]
[[ $(sed -n '3p' "$tmp/report" | cut -f1) == SRR1746971 ]]

source "$repo/scripts/classify_ena_fastq.sh"
while IFS=$'\t' read -r run layout urls md5s; do
  classify_ena_fastq "$run" "$layout" "$urls" "$md5s" 2>"$tmp/classifier.err"
  [[ $FASTQ_LAYOUT == PE && $FASTQ_R1_URL == *"${run}_1.fastq.gz" && $FASTQ_R2_URL == *"${run}_2.fastq.gz" ]]
  [[ $FASTQ_R1_MD5 == m1 && $FASTQ_R2_MD5 == m2 ]]
done < <(tail -n +2 "$tmp/report")

MODE=direct_run CURL_BIN="$mock" "$repo/scripts/resolve_read_runs.sh" SRR1 "$tmp/direct"
[[ $(tail -n +2 "$tmp/direct" | wc -l) -eq 1 ]]
MODE=direct_sample CURL_BIN="$mock" "$repo/scripts/resolve_read_runs.sh" ERS1 "$tmp/sample"
[[ $(tail -n +2 "$tmp/sample" | wc -l) -eq 1 ]]
MODE=mixed_discovery CURL_BIN="$mock" "$repo/scripts/resolve_read_runs.sh" ERS_MIXED "$tmp/mixed"
[[ $(tail -n +2 "$tmp/mixed" | wc -l) -eq 1 ]]
[[ $(sed -n '2p' "$tmp/mixed" | cut -f1) == SRR1 ]]

# ERS sample discovery may return structurally incomplete rows. Both ENA
# discovery layers must remain non-terminal so NCBI can resolve exact runs.
MODE=incomplete_valid CURL_BIN="$mock" "$repo/scripts/resolve_read_runs.sh" ERS14600391 "$tmp/incomplete" 2>"$tmp/incomplete.err"
[[ $(tail -n +2 "$tmp/incomplete" | wc -l) -eq 2 ]]
[[ $(sed -n '2p' "$tmp/incomplete" | cut -f1) == SRR1746970 ]]
grep -q 'no usable exact-run rows; continuing fallback resolution' "$tmp/incomplete.err"

# Once NCBI supplies an exact run, incomplete ENA metadata is terminal.
set +e
MODE=incomplete_strict CURL_BIN="$mock" "$repo/scripts/resolve_read_runs.sh" ERS14600391 "$tmp/strict" 2>"$tmp/strict.err"
strict_rc=$?
set -e
[[ $strict_rc -eq 42 ]]
grep -q 'class=MALFORMED_ENA_METADATA' "$tmp/strict.err"
grep -q 'malformed_exact_run_ena_metadata_SRR999999' "$tmp/strict.err"

set +e
MODE=none CURL_BIN="$mock" "$repo/scripts/resolve_read_runs.sh" SAMN0 "$tmp/none" 2>"$tmp/none.err"
rc=$?
set -e
[[ $rc -eq 42 ]]
grep -q 'class=NO_READ_RUN_RESOLUTION' "$tmp/none.err"
echo 'read-run resolver tests: PASS'
