#!/usr/bin/env bash
# Canonical ENA-to-SAM platform mapping and alignment-header validation.
set -euo pipefail

sam_platform_for_ena() {
    case "$1" in
        ILLUMINA) printf 'ILLUMINA\n' ;;
        DNBSEQ|BGISEQ) printf 'DNBSEQ\n' ;;
        *) echo "ERROR: unsupported ENA instrument platform for SAM read group: $1" >&2; return 2 ;;
    esac
}

derive_sam_platform() {
    [[ $# -eq 1 ]] || { echo "usage: derive_sam_platform BAM_OR_CRAM" >&2; return 2; }
    local input=$1 platforms count platform
    platforms=$(samtools view -H "$input" | awk -F '\t' '
      $1 == "@RG" { for (i=2; i<=NF; i++) if ($i ~ /^PL:/ && length($i)>3) print substr($i,4) }
    ' | LC_ALL=C sort -u)
    count=$(printf '%s\n' "$platforms" | awk 'NF { n++ } END { print n+0 }')
    [[ $count -eq 1 ]] || { echo "ERROR: expected exactly one non-empty @RG PL in $input; found $count: ${platforms:-none}" >&2; return 1; }
    platform=$platforms
    case "$platform" in
        ILLUMINA|DNBSEQ) printf '%s\n' "$platform" ;;
        *) echo "ERROR: unsupported canonical @RG PL in $input: $platform" >&2; return 1 ;;
    esac
}

if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    [[ $# -eq 2 ]] || { echo "usage: $0 map ENA_PLATFORM | derive BAM_OR_CRAM" >&2; exit 2; }
    case "$1" in
        map) sam_platform_for_ena "$2" ;;
        derive) derive_sam_platform "$2" ;;
        *) echo "ERROR: unknown operation: $1" >&2; exit 2 ;;
    esac
fi
