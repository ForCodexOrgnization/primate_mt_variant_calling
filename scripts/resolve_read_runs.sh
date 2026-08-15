#!/usr/bin/env bash
# Resolve an ENA/INSDC alias or NCBI BioSample to ENA read_run metadata.
set -euo pipefail

[[ $# -eq 2 ]] || { echo "usage: $0 INPUT_ACCESSION OUTPUT_TSV" >&2; exit 2; }
accession=$1
output=$2
curl_bin=${CURL_BIN:-curl}
ena_base=${ENA_API_BASE:-https://www.ebi.ac.uk/ena/portal/api}
ncbi_base=${NCBI_EUTILS_BASE:-https://eutils.ncbi.nlm.nih.gov/entrez/eutils}
header=$'run_accession\tlibrary_layout\tfastq_ftp\tfastq_md5'
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

fail_metadata() {
    printf 'ERROR: DETERMINISTIC_FASTQ_FAILURE class=MALFORMED_ENA_METADATA run=%s reason=%s\n' "$accession" "$1" >&2
    exit 42
}
no_resolution() {
    printf 'ERROR: DETERMINISTIC_FASTQ_FAILURE class=NO_READ_RUN_RESOLUTION run=%s reason=no_read_runs_found_after_ena_and_ncbi_queries\n' "$accession" >&2
    exit 42
}
fail_sra_metadata() {
    printf 'ERROR: DETERMINISTIC_FASTQ_FAILURE class=MALFORMED_SRA_METADATA run=%s reason=%s\n' "$accession" "$1" >&2
    exit 42
}
validate_report() {
    local report=$1 context=$2 line run layout urls md5s extra
    [[ -s "$report" ]] || return 1
    [[ "$(head -n1 "$report" | tr -d '\r')" == "$header" ]] || fail_metadata "unexpected_${context}_header"
    while IFS= read -r line; do
        line=${line%$'\r'}
        [[ -n "$line" ]] || continue
        IFS=$'\t' read -r run layout urls md5s extra <<<"$line"
        [[ -z "${extra:-}" && "$run" =~ ^(SRR|ERR|DRR)[0-9]+$ && -n "$layout" && -n "$urls" && -n "$md5s" ]] || fail_metadata "structurally_invalid_${context}_row"
    done < <(tail -n +2 "$report")
    [[ -n "$(tail -n +2 "$report" | sed '/^[[:space:]]*$/d' | head -n1)" ]]
}
ena_report() {
    local acc=$1 dest=$2
    "$curl_bin" -fsSLG "$ena_base/filereport" \
        --data-urlencode "accession=$acc" --data-urlencode result=read_run \
        --data-urlencode fields=run_accession,library_layout,fastq_ftp,fastq_md5 \
        --data-urlencode format=tsv -o "$dest"
}
write_normalized() {
    local source=$1
    { printf '%s\n' "$header"; tail -n +2 "$source" | LC_ALL=C sort -t$'\t' -k1,1 -k2,2 | awk -F '\t' '!seen[$1]++'; } >"$output"
}

# A transport failure is transient (non-42).  A valid header with no rows is a
# normal miss and advances to the next resolver layer.
set +e; ena_report "$accession" "$tmp/direct.tsv"; direct_rc=$?; set -e
[[ $direct_rc -eq 0 || $direct_rc -eq 22 ]] || exit "$direct_rc"
if [[ $direct_rc -eq 0 ]] && validate_report "$tmp/direct.tsv" direct_ena; then
    write_normalized "$tmp/direct.tsv"
    exit 0
fi

query="run_accession=\"$accession\" OR experiment_accession=\"$accession\" OR study_accession=\"$accession\" OR secondary_study_accession=\"$accession\" OR sample_accession=\"$accession\" OR secondary_sample_accession=\"$accession\" OR sample_alias=\"$accession\""
set +e
"$curl_bin" -fsSLG "$ena_base/search" --data-urlencode result=read_run \
    --data-urlencode "query=$query" \
    --data-urlencode fields=run_accession,library_layout,fastq_ftp,fastq_md5 \
    --data-urlencode format=tsv -o "$tmp/search.tsv"
search_rc=$?
set -e
[[ $search_rc -eq 0 || $search_rc -eq 22 ]] || exit "$search_rc"
if [[ $search_rc -eq 0 ]] && validate_report "$tmp/search.tsv" advanced_ena; then
    write_normalized "$tmp/search.tsv"
    exit 0
fi

"$curl_bin" -fsSLG "$ncbi_base/esearch.fcgi" --data-urlencode db=sra \
    --data-urlencode "term=$accession" --data-urlencode retmode=json -o "$tmp/esearch.json"
set +e
python3 - "$tmp/esearch.json" "$tmp/ids" <<'PY'
import json, sys
try:
    data = json.load(open(sys.argv[1]))
    ids = data["esearchresult"]["idlist"]
    if not isinstance(ids, list): raise ValueError
except (KeyError, TypeError, ValueError, json.JSONDecodeError):
    sys.exit(42)
open(sys.argv[2], "w").write(",".join(str(x) for x in ids))
PY
json_rc=$?
set -e
case $json_rc in 0) :;; 42) fail_sra_metadata malformed_ncbi_esearch_response;; *) exit "$json_rc";; esac
[[ -s "$tmp/ids" ]] || no_resolution
"$curl_bin" -fsSLG "$ncbi_base/efetch.fcgi" --data-urlencode db=sra \
    --data-urlencode "id=$(cat "$tmp/ids")" --data-urlencode retmode=xml -o "$tmp/sra.xml"
set +e
python3 - "$tmp/sra.xml" "$tmp/runs" <<'PY'
import re, sys, xml.etree.ElementTree as ET
from pathlib import Path
text = Path(sys.argv[1]).read_text(errors="strict")
root = ET.fromstring(text)
runs = sorted({value for elem in root.iter() for value in elem.attrib.values()
               if re.fullmatch(r'(?:SRR|ERR|DRR)\d+', value)})
Path(sys.argv[2]).write_text("\n".join(runs) + ("\n" if runs else ""))
PY
xml_rc=$?
set -e
[[ $xml_rc -eq 0 ]] || fail_sra_metadata malformed_ncbi_sra_xml
[[ -s "$tmp/runs" ]] || no_resolution

printf '%s\n' "$header" >"$output"
while IFS= read -r run; do
    ena_report "$run" "$tmp/$run.tsv"
    validate_report "$tmp/$run.tsv" resolved_run_ena || fail_metadata "resolved_run_missing_ena_metadata_${run}"
    [[ $(tail -n +2 "$tmp/$run.tsv" | sed '/^[[:space:]]*$/d' | wc -l) -eq 1 ]] || fail_metadata "resolved_run_nonunique_ena_metadata_${run}"
    [[ $(tail -n +2 "$tmp/$run.tsv" | cut -f1) == "$run" ]] || fail_metadata "resolved_run_accession_mismatch_${run}"
    tail -n +2 "$tmp/$run.tsv" >>"$output"
done <"$tmp/runs"

# NCBI packages can mention the same run in several XML attributes; enforce a
# final deterministic run order without touching the semicolon-delimited URL
# and MD5 vectors inside each row.
{ printf '%s\n' "$header"; tail -n +2 "$output" | LC_ALL=C sort -t$'\t' -k1,1; } >"$tmp/final.tsv"
mv "$tmp/final.tsv" "$output"
