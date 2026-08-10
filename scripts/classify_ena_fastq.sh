#!/usr/bin/env bash

# Classify one ENA read_run row without relying on the order of its FASTQ URLs.
# On success this function sets FASTQ_* globals.  On deterministic metadata or
# layout errors it prints a machine-readable diagnostic and returns 42.
classify_ena_fastq() {
    local run_id=$1 metadata_layout=$2 ftp_urls=$3 md5s=$4
    local -a urls mds
    local i name r1_count=0 r2_count=0 orphan_count=0 other_count=0

    FASTQ_LAYOUT= FASTQ_DETECTED_LAYOUT= FASTQ_R1_URL= FASTQ_R1_MD5=
    FASTQ_R2_URL= FASTQ_R2_MD5= FASTQ_ORPHAN_URL= FASTQ_ORPHAN_MD5=

    deterministic_fastq_error() {
        local class=$1 reason=$2
        printf 'ERROR: DETERMINISTIC_FASTQ_FAILURE class=%s run=%s reason=%s\n' \
            "$class" "$run_id" "$reason" >&2
        return 42
    }

    [[ -n "$run_id" && -n "$metadata_layout" && -n "$ftp_urls" && -n "$md5s" ]] || {
        deterministic_fastq_error MALFORMED_ENA_METADATA "missing required ENA read_run field"
        return 42
    }
    IFS=';' read -r -a urls <<< "$ftp_urls"
    IFS=';' read -r -a mds <<< "$md5s"
    if (( ${#urls[@]} != ${#mds[@]} )); then
        deterministic_fastq_error MD5_URL_COUNT_MISMATCH "url_count=${#urls[@]} md5_count=${#mds[@]}"
        return 42
    fi
    (( ${#urls[@]} > 0 )) || { deterministic_fastq_error MALFORMED_ENA_METADATA "no FASTQ URLs"; return 42; }

    # URL and checksum are inspected and selected at the same index.  They are
    # deliberately never sorted into independent arrays.
    for i in "${!urls[@]}"; do
        [[ -n "${urls[$i]}" && -n "${mds[$i]}" ]] || {
            deterministic_fastq_error MALFORMED_ENA_METADATA "empty URL or MD5 at index $i"
            return 42
        }
        name=$(basename "${urls[$i]}")
        case "$name" in
            "${run_id}_1.fastq.gz")
                ((++r1_count)); FASTQ_R1_URL=${urls[$i]}; FASTQ_R1_MD5=${mds[$i]} ;;
            "${run_id}_2.fastq.gz")
                ((++r2_count)); FASTQ_R2_URL=${urls[$i]}; FASTQ_R2_MD5=${mds[$i]} ;;
            "${run_id}.fastq.gz")
                ((++orphan_count)); FASTQ_ORPHAN_URL=${urls[$i]}; FASTQ_ORPHAN_MD5=${mds[$i]} ;;
            *) ((++other_count)) ;;
        esac
    done

    (( r1_count <= 1 )) || { deterministic_fastq_error DUPLICATE_R1 "multiple ${run_id}_1.fastq.gz entries"; return 42; }
    (( r2_count <= 1 )) || { deterministic_fastq_error DUPLICATE_R2 "multiple ${run_id}_2.fastq.gz entries"; return 42; }
    (( orphan_count <= 1 )) || { deterministic_fastq_error AMBIGUOUS_FASTQ_LAYOUT "duplicate ${run_id}.fastq.gz entries"; return 42; }
    (( other_count == 0 )) || { deterministic_fastq_error AMBIGUOUS_FASTQ_LAYOUT "unrecognized FASTQ basename count=$other_count"; return 42; }

    if (( r1_count == 1 && r2_count == 1 )); then
        if (( orphan_count == 1 && ${#urls[@]} == 3 )); then
            FASTQ_DETECTED_LAYOUT=PE_PLUS_ORPHAN
            echo "WARN: run=$run_id archive layout contains paired R1/R2 plus orphan FASTQ; using R1/R2 and ignoring orphan FASTQ." >&2
            if [[ "$metadata_layout" != PAIRED ]]; then
                echo "WARN: ENA_METADATA_LAYOUT_MISMATCH run=$run_id metadata_layout=$metadata_layout detected_layout=PE_PLUS_ORPHAN action=USE_R1_R2_IGNORE_ORPHAN" >&2
            fi
        elif (( orphan_count == 0 && ${#urls[@]} == 2 )); then
            FASTQ_DETECTED_LAYOUT=PE
            [[ "$metadata_layout" == PAIRED ]] || echo "WARN: ENA_METADATA_LAYOUT_MISMATCH run=$run_id metadata_layout=$metadata_layout detected_layout=PE action=USE_R1_R2" >&2
        else
            deterministic_fastq_error AMBIGUOUS_FASTQ_LAYOUT "paired FASTQs accompanied by unexplained entries"
            return 42
        fi
        FASTQ_LAYOUT=PE
        return 0
    fi

    if (( r1_count == 1 && r2_count == 0 && orphan_count == 0 && ${#urls[@]} == 1 )); then
        FASTQ_LAYOUT=SE; FASTQ_DETECTED_LAYOUT=SE_UNUSUAL_R1_NAME
        echo "WARN: run=$run_id single FASTQ uses unusual ${run_id}_1.fastq.gz filename; treating as SE." >&2
        return 0
    fi
    (( r1_count == 0 && r2_count == 0 && orphan_count == 1 && ${#urls[@]} == 1 )) && {
        FASTQ_LAYOUT=SE; FASTQ_DETECTED_LAYOUT=SE
        FASTQ_R1_URL=$FASTQ_ORPHAN_URL; FASTQ_R1_MD5=$FASTQ_ORPHAN_MD5
        return 0
    }
    (( r1_count == 0 && r2_count == 1 )) && { deterministic_fastq_error MISSING_R1 "R2 exists without R1"; return 42; }
    (( r1_count == 1 && r2_count == 0 )) && { deterministic_fastq_error MISSING_R2 "R1 exists without R2 among multiple FASTQs"; return 42; }
    deterministic_fastq_error UNSUPPORTED_FASTQ_LAYOUT "FASTQ layout could not be classified"
    return 42
}
