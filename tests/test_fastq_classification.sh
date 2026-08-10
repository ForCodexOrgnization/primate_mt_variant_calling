#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$repo/scripts/classify_ena_fastq.sh"

ok() {
    local run=$1 metadata=$2 urls=$3 md5s=$4 expected_layout=$5 expected_detected=$6
    classify_ena_fastq "$run" "$metadata" "$urls" "$md5s" 2>"$tmp/stderr"
    [[ "$FASTQ_LAYOUT" == "$expected_layout" && "$FASTQ_DETECTED_LAYOUT" == "$expected_detected" ]]
}
fails_as() {
    local expected=$1 run=$2 metadata=$3 urls=$4 md5s=$5
    if classify_ena_fastq "$run" "$metadata" "$urls" "$md5s" 2>"$tmp/stderr"; then
        echo "classification unexpectedly succeeded: $expected" >&2; exit 1
    else
        rc=$?
    fi
    [[ $rc == 42 ]]
    grep -F "class=$expected" "$tmp/stderr" >/dev/null
}

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

ok ERR1 SINGLE 'host/ERR1.fastq.gz' a SE SE
[[ "$FASTQ_R1_URL:$FASTQ_R1_MD5" == 'host/ERR1.fastq.gz:a' ]]

# Reverse order proves R1/R2 identity is filename-derived, not array-derived.
ok ERR2 PAIRED 'host/ERR2_2.fastq.gz;host/ERR2_1.fastq.gz' 'b;a' PE PE
[[ "$FASTQ_R1_URL:$FASTQ_R1_MD5" == 'host/ERR2_1.fastq.gz:a' ]]
[[ "$FASTQ_R2_URL:$FASTQ_R2_MD5" == 'host/ERR2_2.fastq.gz:b' ]]

ok ERR3 PAIRED 'host/ERR3.fastq.gz;host/ERR3_2.fastq.gz;host/ERR3_1.fastq.gz' 'o;b;a' PE PE_PLUS_ORPHAN
grep -F 'ignoring orphan FASTQ' "$tmp/stderr" >/dev/null
[[ "$FASTQ_ORPHAN_URL" == host/ERR3.fastq.gz ]]

# Real regression fixture: ENA says SINGLE, archive filenames contain a pair
# plus an orphan.  No actual archive data is downloaded.
ok ERR11028341 SINGLE \
  'ftp.sra.ebi.ac.uk/vol1/fastq/ERR110/041/ERR11028341/ERR11028341.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR110/041/ERR11028341/ERR11028341_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR110/041/ERR11028341/ERR11028341_2.fastq.gz' \
  'a3720b070746b13ef5da16cb105d844a;c7f0d9947267cdc8b594fa6021021d89;558eac0803d77bfc0835dd8be1d281af' PE PE_PLUS_ORPHAN
grep -F 'ENA_METADATA_LAYOUT_MISMATCH run=ERR11028341 metadata_layout=SINGLE detected_layout=PE_PLUS_ORPHAN action=USE_R1_R2_IGNORE_ORPHAN' "$tmp/stderr" >/dev/null
[[ "$FASTQ_R1_MD5" == c7f0d9947267cdc8b594fa6021021d89 ]]
[[ "$FASTQ_R2_MD5" == 558eac0803d77bfc0835dd8be1d281af ]]

fails_as MD5_URL_COUNT_MISMATCH ERR4 SINGLE 'h/ERR4.fastq.gz' 'a;b'
fails_as DUPLICATE_R1 ERR5 PAIRED 'a/ERR5_1.fastq.gz;b/ERR5_1.fastq.gz;c/ERR5_2.fastq.gz' 'a;b;c'
fails_as DUPLICATE_R2 ERR6 PAIRED 'a/ERR6_1.fastq.gz;b/ERR6_2.fastq.gz;c/ERR6_2.fastq.gz' 'a;b;c'
fails_as AMBIGUOUS_FASTQ_LAYOUT ERR7 PAIRED 'a/ERR7.fastq.gz;a/ERR7_1.fastq.gz;a/ERR7_2.fastq.gz;a/unrelated.fastq.gz' 'a;b;c;d'

echo 'FASTQ classification tests: PASS'
