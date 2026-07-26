#!/usr/bin/env nextflow
nextflow.enable.dsl=2

if (!params.containsKey('skip_existing_cram') || params.skip_existing_cram == null) {
    params.skip_existing_cram = true
}
params.force_reprocess_existing_cram = params.force_reprocess_existing_cram ?: false
params.backfill_existing_cram_marker = params.containsKey('backfill_existing_cram_marker') ? params.backfill_existing_cram_marker : true
params.existing_cram_check_retries = params.existing_cram_check_retries ?: 3
params.existing_cram_check_delay_seconds = params.existing_cram_check_delay_seconds ?: 10
params.existing_cram_min_size_bytes = params.existing_cram_min_size_bytes ?: 1024
params.existing_crai_min_size_bytes = params.existing_crai_min_size_bytes ?: 16
params.existing_cram_quickcheck_retries = params.existing_cram_quickcheck_retries ?: 3
params.existing_cram_quickcheck_delay_seconds = params.existing_cram_quickcheck_delay_seconds ?: 10
params.existing_cram_quickcheck_timeout_seconds = params.existing_cram_quickcheck_timeout_seconds ?: 120
def paramAsBoolean = { value -> value instanceof Boolean ? value : value?.toString()?.toBoolean() }

if (!params.containsKey('enable_chunked_alignment') || params.enable_chunked_alignment == null) {
    params.enable_chunked_alignment = true
}
if (!params.containsKey('chunked_alignment_size_threshold_gb') || params.chunked_alignment_size_threshold_gb == null) {
    params.chunked_alignment_size_threshold_gb = 20
}
if (!params.containsKey('fastq_chunk_reads') || params.fastq_chunk_reads == null) {
    params.fastq_chunk_reads = 25000000
}
if (!params.containsKey('cleanup_chunk_intermediates') || params.cleanup_chunk_intermediates == null) {
    params.cleanup_chunk_intermediates = false
}
if (!params.containsKey('align_sort_threads') || params.align_sort_threads == null) {
    params.align_sort_threads = 4
}
if (!params.containsKey('align_sort_mem_per_thread') || params.align_sort_mem_per_thread == null) {
    params.align_sort_mem_per_thread = '2G'
}

def runCommand = { List command ->
    def stringCommand = command.collect { item ->
        if (item == null) {
            throw new IllegalArgumentException("Null command argument in: ${command}")
        }
        item.toString()
    }

    try {
        log.info "Running external command: ${stringCommand.collect { it.contains(' ') ? "'${it}'" : it }.join(' ')}"

        def proc = new ProcessBuilder(stringCommand)
            .redirectErrorStream(true)
            .start()
        def output = new StringBuffer()
        def reader = Thread.start {
            proc.inputStream.withReader { stream ->
                stream.eachLine { line ->
                    output.append(line).append('\n')
                }
            }
        }

        proc.waitFor()
        reader.join(5000)

        return [exitCode: proc.exitValue(), output: output.toString(), exception: null]
    } catch (Exception e) {
        return [exitCode: null, output: "${e.class.name}: ${e.message}", exception: e]
    }
}

def validateExistingCram = { String sampleId, String refName ->
    def alignmentDir = new File("${params.outdir}/${sampleId}/alignment")
    def cram = new File(alignmentDir, "${sampleId}.cram")
    def candidates = [new File(alignmentDir, "${sampleId}.cram.crai"), new File(alignmentDir, "${sampleId}.crai")]
    def crai = candidates.find { it.isFile() }
    if (!cram.isFile()) return [status:'INCOMPLETE', reason:'cram_missing', cram:cram.absolutePath, crai:null]
    if (!crai) return [status:'INCOMPLETE', reason:"crai_missing_checked=${candidates*.absolutePath.join(',')}", cram:cram.absolutePath, crai:null]

    def marker = new File(alignmentDir, "${sampleId}.cram.complete")
    def ref = new File("${params.global_ref_dir}/${refName}.fa")
    def validatorScript = new File(projectDir.toString(), 'scripts/validate_cram.sh').absolutePath
    def command = [validatorScript, '--cram', cram.absolutePath, '--crai', crai.absolutePath,
        '--reference', ref.absolutePath, '--samtools', params.samtools_bin.toString(),
        '--stability-retries', params.existing_cram_check_retries.toString(),
        '--retries', params.existing_cram_quickcheck_retries.toString(),
        '--delay', params.existing_cram_quickcheck_delay_seconds.toString(),
        '--timeout', params.existing_cram_quickcheck_timeout_seconds.toString(),
        '--min-cram-size', params.existing_cram_min_size_bytes.toString(), '--min-crai-size', params.existing_crai_min_size_bytes.toString()
    ].collect { it.toString() }
    if (marker.isFile()) command.addAll(['--marker', marker.absolutePath].collect { it.toString() })
    assert command.every { it instanceof String }
    def result = runCommand(command)
    def status = result.exitCode == 0 ? 'COMPLETE' : (result.exitCode == 1 ? 'INCOMPLETE' : 'UNKNOWN')
    def reasonMatch = (result.output =~ /(?m)^REASON=(.*)$/)
    def reason = reasonMatch.find() ? reasonMatch.group(1) : "validator_exit=${result.exitCode}"
    def validation = [status:status, reason:reason, cram:cram.absolutePath, crai:crai.absolutePath,
        marker:marker.isFile() ? 'present' : 'legacy_missing', diagnostics:result.output.trim()]

    if (status == 'COMPLETE' && !marker.isFile() && paramAsBoolean(params.backfill_existing_cram_marker)) {
        try {
            def tmp = new File(alignmentDir, ".${sampleId}.cram.complete.${UUID.randomUUID()}.tmp")
            tmp.text = "sample_id=${sampleId}\nref_name=${refName}\ncram=${cram.name}\ncrai=${crai.name}\ncram_size=${cram.length()}\ncrai_size=${crai.length()}\nsamtools_version=coordinator_validation\ncompleted_at=${java.time.OffsetDateTime.now()}\n"
            java.nio.file.Files.move(tmp.toPath(), marker.toPath(), java.nio.file.StandardCopyOption.ATOMIC_MOVE)
            validation.marker = 'backfilled'
        } catch (Exception e) {
            log.warn "Validated legacy CRAM for ${sampleId}, but atomic marker backfill failed: ${e}"
        }
    }
    validation
}

/*
================================================================================
    STRICT PRIMATE ALIGNMENT PIPELINE
================================================================================
*/

log.info """
STRICT PRIMATE ALIGNMENT PIPELINE START
======================================
Sample TSV          : ${params.sample_tsv}
Output Directory    : ${params.outdir}
Reference Directory : ${params.global_ref_dir}
Skip existing CRAM  : ${params.skip_existing_cram}
Force reprocess     : ${params.force_reprocess_existing_cram}
Samtools binary     : ${params.samtools_bin}
Chunked alignment   : ${params.enable_chunked_alignment}
Chunk threshold GB  : ${params.chunked_alignment_size_threshold_gb}
FASTQ chunk reads   : ${params.fastq_chunk_reads}
======================================
"""

// 1. 从 TSV 读取初始样本信息
//
// Important: Nextflow prints process names with "[-]" when their input channels are
// empty.  That looks like a skip even when no complete CRAM was found.  Fail fast
// for empty/malformed sample batches so an empty output directory cannot be
// reported as a successful run.
ch_parsed_samples = Channel.fromPath(params.sample_tsv)
    .ifEmpty { error "Sample TSV path did not match any file: ${params.sample_tsv}" }
    .splitCsv(header: false, sep: '\t', strip: true)
    .filter { row ->
        if (row.size() < 2 || !row[0]?.trim()) {
            log.warn "Ignoring malformed/blank sample TSV row: ${row}"
            return false
        }
        return true
    }
    .filter { row ->
        def sample_id = row[0].trim()
        if (sample_id.equalsIgnoreCase('sample') || sample_id.equalsIgnoreCase('sample_id')) {
            log.warn "Ignoring sample TSV header row: ${row}"
            return false
        }
        return true
    }
    .map { row ->
        def meta = [id: row[0].trim()]
        def species = row[1].trim()
        // Accept both legacy 2-column batches (sample_id, species_name) and
        // 3-column batches (sample_id, species_name, ref_name).  For the
        // current reference layout, ref_name defaults to species_name.
        def ref_name = row.size() >= 3 && row[2]?.trim() ? row[2].trim() : species
        tuple(meta, species, ref_name)
    }
    .ifEmpty { error "No valid samples found in ${params.sample_tsv}; expected tab-separated rows: sample_id<TAB>species_name[<TAB>ref_name]" }

ch_samples = ch_parsed_samples
    // Preserve each emitted tuple as one list element. Nextflow's collect()
    // flattens tuple items by default, which would turn tuples into a stream
    // like [meta, species, ref_name, ...] and make sample_tuple[0] null for
    // the meta map entries below.
    .collect(flat: false)
    .flatMap { parsed_samples ->
        def samples_to_process = parsed_samples.findAll { sample_tuple ->
            if (!sample_tuple || sample_tuple.size() < 3 || !sample_tuple[0]?.id) {
                error "Malformed parsed sample entry in ${params.sample_tsv}: ${sample_tuple}"
            }

            def meta = sample_tuple[0]
            if (paramAsBoolean(params.force_reprocess_existing_cram)) {
                log.warn "EXISTING_CRAM_CHECK sample=${meta.id} action=reprocess reason=user_forced"
                return true
            }
            if (!paramAsBoolean(params.skip_existing_cram)) return true
            def validation = validateExistingCram(meta.id, sample_tuple[2])
            def action = validation.status == 'COMPLETE' ? 'skip' : (validation.status == 'INCOMPLETE' ? 'reprocess' : 'abort_not_reprocess')
            log.info "EXISTING_CRAM_CHECK sample=${meta.id} status=${validation.status} cram=${validation.cram} crai=${validation.crai} marker=${validation.marker ?: 'absent'} reason=${validation.reason} action=${action}\n${validation.diagnostics ?: ''}"
            if (validation.status == 'COMPLETE') return false
            if (validation.status == 'INCOMPLETE') return true
            if (validation.status == 'UNKNOWN') error "Cannot determine whether existing CRAM for sample ${meta.id} is valid. ${validation.reason}. Refusing to re-download and realign; resolve coordinator/storage validation or use --force_reprocess_existing_cram true."
            error "Unexpected existing CRAM validation status: ${validation}"
            return true
        }

        if (samples_to_process.isEmpty()) {
            log.info "All samples in this batch already have complete CRAM/CRAI; nothing to run in preprocessing."
        }

        return samples_to_process
    }

workflow {

    // 第一步：下载 FASTQ 并生成配对清单
    DOWNLOAD_FASTQ(ch_samples)

    // 第二步：解析下载生成的 fastq_pairs.tsv，展平为单对任务流
    ch_fastq_pairs = DOWNLOAD_FASTQ.out.pairs_tsv
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [
                id      : row.sample_id,
                pair_id : row.run_id,
                layout  : row.layout,
                n_pairs : row.expected_pairs.toInteger()
            ]
            // 必须使用 file() 包装路径，Nextflow 才能在进程间正确传递文件
            def reads = row.layout == 'PE' ? [ file(row.r1), file(row.r2) ] : [ file(row.r1) ]
            tuple(meta, row.species_name, row.ref_name, reads)
        }

    // 第三步：根据 FASTQ 总大小路由到常规比对或分块比对
    ch_fastq_alignment_routes = ch_fastq_pairs
        .map { meta, species_name, ref_name, reads ->
            long total_fastq_bytes = reads.collect { it.size() }.sum() as long
            double total_fastq_size_gb = total_fastq_bytes / 1024.0 / 1024.0 / 1024.0
            boolean use_chunked = params.enable_chunked_alignment && total_fastq_size_gb >= (params.chunked_alignment_size_threshold_gb as double)
            String route_reason = !params.enable_chunked_alignment ? 'user_disabled' : (use_chunked ? 'above_threshold' : 'below_threshold')
            def routed_meta = meta + [
                total_fastq_size_gb: total_fastq_size_gb,
                use_chunked_alignment: use_chunked
            ]
            log.info String.format(
                "FASTQ sizing: sample=%s run=%s layout=%s total_fastq_size_gb=%.3f threshold_gb=%s chunked_enabled=%s route=%s reason=%s",
                meta.id,
                meta.pair_id,
                meta.layout,
                total_fastq_size_gb,
                params.chunked_alignment_size_threshold_gb,
                params.enable_chunked_alignment,
                use_chunked ? 'chunked' : 'standard',
                route_reason
            )
            tuple(routed_meta, species_name, ref_name, reads)
        }
        .branch { meta, species_name, ref_name, reads ->
            chunked: meta.use_chunked_alignment
            standard: true
        }

    ALIGN_AND_SORT(ch_fastq_alignment_routes.standard)

    // Make each large-run chunk an independently cached Nextflow task. The old
    // streaming process remains below only as a deprecated manual fallback.
    SPLIT_FASTQ_CHUNKS(ch_fastq_alignment_routes.chunked)

    ch_alignment_chunks = SPLIT_FASTQ_CHUNKS.out.chunks_tsv
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [
                id       : row.sample_id,
                pair_id  : row.run_id,
                layout   : row.layout,
                chunk_id : row.chunk_id,
                n_chunks : row.expected_chunks.toInteger(),
                n_pairs  : row.expected_pairs.toInteger()
            ]
            def reads = row.layout == 'PE' ? [file(row.r1), file(row.r2)] : [file(row.r1)]
            tuple(meta, row.species_name, row.ref_name, reads)
        }

    ALIGN_AND_SORT_CHUNK(ch_alignment_chunks)

    ch_run_chunk_bams = ALIGN_AND_SORT_CHUNK.out.bam
        .map { meta, species_name, ref_name, bam, bai ->
            tuple(groupKey([sample_id: meta.id, pair_id: meta.pair_id], meta.n_chunks),
                  meta, species_name, ref_name, bam)
        }
        .groupTuple()
        .map { run_key, metas, species_names, ref_names, bams ->
            def meta = metas[0]
            if (bams.size() != meta.n_chunks) {
                error "Missing chunk BAMs for ${meta.id}.${meta.pair_id}: expected ${meta.n_chunks}, received ${bams.size()}"
            }
            if (metas.any { it.n_chunks != meta.n_chunks || it.n_pairs != meta.n_pairs }) {
                error "Inconsistent chunk metadata for ${meta.id}.${meta.pair_id}"
            }
            tuple(meta, species_names[0], ref_names[0], bams.sort { it.name })
        }

    MERGE_RUN_CHUNK_BAMS(ch_run_chunk_bams)

    ch_per_run_bams = ALIGN_AND_SORT.out.bam.mix(MERGE_RUN_CHUNK_BAMS.out.bam)

    // 第四步：按 Sample ID 分组，并强制校验数量
    ch_bam_grouped = ch_per_run_bams
        .map { meta, species, ref_name, bam, bai ->
            // groupKey 确保收齐 n_pairs 个文件后才下发到下游
            tuple(groupKey(meta.id, meta.n_pairs), meta, species, ref_name, bam)
        }
        .groupTuple()
        .map { gKey, metas, species_list, refs, bams ->
            def meta = metas[0]
            // 严苛性检查：如果实际 BAM 数量不等于预期数量，报错终止
            if (bams.size() != meta.n_pairs) {
                error "CRITICAL: Sample ${meta.id} missing runs. Expected ${meta.n_pairs}, got ${bams.size()}."
            }
            tuple(meta, species_list[0], refs[0], bams)
        }

    // 第五步：合并样本的所有 BAM 文件
    MERGE_BAMS(ch_bam_grouped)

    // 第六步：转为 CRAM 并发布结果
    BAM_TO_CRAM(MERGE_BAMS.out.merged_bam)
}

/*
================================================================================
    PROCESSES
================================================================================
*/

process DOWNLOAD_FASTQ {
    tag "${meta.id}"
    label 'down_task'
    
    // Retry transient task-level failures, then skip only the bad sample so the batch can continue.
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(meta), val(species_name), val(ref_name)

    output:
    path "fastq_pairs.tsv", emit: pairs_tsv
    path "failed_download.tsv", emit: failed_tsv, optional: true
    path "fastqs/*.fastq.gz", emit: fastq_files

    script:
    // 将你的 aria2c 路径定义为 Nextflow 变量
    def aria2_bin = "/home/lt692/.conda/envs/aria2_env/bin/aria2c"
    def aria2_connections = params.aria2_connections as int
    def aria2_splits = params.aria2_splits as int
    def aria2_fallback_connections = params.aria2_fallback_connections as int
    def aria2_fallback_splits = params.aria2_fallback_splits as int
    def aria2_fallback_retry_wait = params.aria2_fallback_retry_wait as int
    def max_download_attempts = params.max_download_attempts as int
    
    """
    #!/usr/bin/env bash
    set -Eeuo pipefail

    cleanup_download_outputs() {
        echo "WARN: DOWNLOAD_FASTQ interrupted or failed for ${meta.id}; keeping downloaded FASTQ files for possible resume." >&2
        rm -f report.tsv fastq_pairs.tsv failed_download.tsv
    }

    trap cleanup_download_outputs INT TERM HUP QUIT ERR

    mkdir -p fastqs

    # Do not remove fastqs here. Keep completed or partial downloads so aria2c -c can resume.
    rm -f report.tsv fastq_pairs.tsv failed_download.tsv
    ena_base="https://www.ebi.ac.uk/ena/portal/api"
    acc="${meta.id}"

    # 1. 获取 Run IDs
    # ENA filereport accepts INSDC accessions, but some sample sheets use aliases
    # such as HG00119. If filereport returns HTTP 400, fall back to the Portal
    # advanced search API and resolve the alias to read_run records.
    if ! curl -fsSLG "\${ena_base}/filereport" \
        --data-urlencode "accession=\${acc}" \
        --data-urlencode "result=read_run" \
        --data-urlencode "fields=run_accession,fastq_ftp,fastq_md5" \
        --data-urlencode "format=tsv" \
        -o report.tsv; then
        echo "WARN: ENA filereport failed for \${acc}; trying read_run search by accession/sample alias." >&2
        ena_query="run_accession=\"\${acc}\" OR experiment_accession=\"\${acc}\" OR study_accession=\"\${acc}\" OR secondary_study_accession=\"\${acc}\" OR sample_accession=\"\${acc}\" OR secondary_sample_accession=\"\${acc}\" OR sample_alias=\"\${acc}\""
        curl -fsSLG "\${ena_base}/search" \
            --data-urlencode "result=read_run" \
            --data-urlencode "query=\${ena_query}" \
            --data-urlencode "fields=run_accession,fastq_ftp,fastq_md5" \
            --data-urlencode "format=tsv" \
            -o report.tsv
    fi

    if [[ \$(tail -n +2 report.tsv | wc -l) -eq 0 ]]; then
        echo "ERROR: No runs found for \${acc}. Checked ENA accession fields and sample_alias." >&2; exit 1
    fi

    # 2. 准备输出清单
    rm -f fastq_pairs.tmp

    record_failed_download() {
        local failed_run_id="\$1"
        local failed_fastq="\$2"
        local reason="\$3"
        if [[ ! -f failed_download.tsv ]]; then
            echo -e "sample_id\trun_id\tfastq_file\treason" > failed_download.tsv
        fi
        echo -e "${meta.id}\t\${failed_run_id}\t\${failed_fastq}\t\${reason}" >> failed_download.tsv
    }

    # 3. 循环下载并验证 MD5
    tail -n +2 report.tsv | while IFS=\$'\\t' read -r run_id ftp_urls md5s; do
        IFS=';' read -r -a urls <<< "\$ftp_urls"
        IFS=';' read -r -a mds <<< "\$md5s"

        if [[ \${#urls[@]} -eq 2 ]]; then
            layout="PE"
            echo "INFO: Run \$run_id detected as paired-end."
        elif [[ \${#urls[@]} -eq 1 ]]; then
            layout="SE"
            echo "INFO: Run \$run_id detected as single-end."
        else
            echo "ERROR: Run \$run_id has unsupported number of FASTQ URLs: \${#urls[@]}." >&2
            record_failed_download "\$run_id" "." "UNSUPPORTED_FASTQ_URL_COUNT_\${#urls[@]}"
            exit 1
        fi

        if [[ \${#mds[@]} -ne \${#urls[@]} ]]; then
            echo "ERROR: Run \$run_id has \${#urls[@]} FASTQ URLs but \${#mds[@]} MD5 values." >&2
            record_failed_download "\$run_id" "." "MD5_URL_COUNT_MISMATCH"
            exit 1
        fi

        download_and_verify_fastq() {
            local url="\$1"
            local target="\$2"
            local expected_md5="\$3"
            local run_id="\$4"
            local fastq_file
            fastq_file=\$(basename "\$target")

            if [[ -z "\$expected_md5" ]]; then
                echo "ERROR: Missing MD5 checksum for \$target in ENA report." >&2
                rm -f "\$target" "\$target.aria2"
                record_failed_download "\$run_id" "\$fastq_file" "MISSING_MD5"
                return 1
            fi

            if [[ -s "\$target" ]]; then
                echo "INFO: Existing FASTQ found for \$target; checking gzip and MD5..."
                if gzip -t "\$target" && echo "\$expected_md5  \$target" | md5sum -c -; then
                    echo "INFO: Existing FASTQ passed gzip and MD5; skipping download for \$target."
                    rm -f "\$target.aria2"
                    return 0
                else
                    echo "WARN: Existing FASTQ failed gzip or MD5 validation; removing and re-downloading \$target." >&2
                    rm -f "\$target" "\$target.aria2"
                fi
            fi

            local aria2_log="\${target}.aria2.log"
            local attempt=1
            while (( attempt <= ${max_download_attempts} )); do
                echo "INFO: Normal aria2c attempt \$attempt/${max_download_attempts} for \$target with -x ${aria2_connections} -s ${aria2_splits}. Log: \$aria2_log"
                if ${aria2_bin} -x ${aria2_connections} -s ${aria2_splits} -c -m 5 --retry-wait 10 \
                    --file-allocation=none \
                    --summary-interval=0 -d fastqs -o "\$fastq_file" "\$url" >"\$aria2_log" 2>&1; then
                    echo "INFO: Normal aria2c attempt completed for \$target; validating gzip and MD5."
                elif grep -Eq 'status=403|status=429|errorCode=22' "\$aria2_log"; then
                    echo "WARN: Normal aria2c attempt for \$target hit HTTP 403/429 or aria2 errorCode=22; retrying in fallback single-connection mode. See \$aria2_log" >&2
                    if ${aria2_bin} -x ${aria2_fallback_connections} -s ${aria2_fallback_splits} -c -m 5 \
                        --retry-wait=${aria2_fallback_retry_wait} \
                        --timeout=120 \
                        --connect-timeout=60 \
                        --file-allocation=none \
                        --summary-interval=0 -d fastqs -o "\$fastq_file" "\$url" >>"\$aria2_log" 2>&1; then
                        echo "INFO: Fallback aria2c single-connection mode completed for \$target; validating gzip and MD5."
                    else
                        echo "WARN: Fallback aria2c single-connection mode failed for \$target. See \$aria2_log" >&2
                        rm -f "\$target" "\$target.aria2"
                        ((attempt++))
                        continue
                    fi
                else
                    echo "WARN: Normal aria2c attempt failed for \$target without a fallback-triggering status. See \$aria2_log" >&2
                    rm -f "\$target" "\$target.aria2"
                    ((attempt++))
                    continue
                fi

                if ! gzip -t "\$target"; then
                    echo "WARN: Final gzip status for \$target: FAILED. Removing corrupted file before retry." >&2
                    rm -f "\$target" "\$target.aria2"
                    ((attempt++))
                    continue
                fi
                echo "INFO: Final gzip status for \$target: PASSED."

                if ! echo "\$expected_md5  \$target" | md5sum -c -; then
                    echo "WARN: Final MD5 status for \$target: FAILED. Removing corrupted file before retry." >&2
                    rm -f "\$target" "\$target.aria2"
                    ((attempt++))
                    continue
                fi

                echo "INFO: Final MD5 status for \$target: PASSED."
                echo "INFO: gzip and MD5 checks passed for \$target."
                rm -f "\$target.aria2"
                return 0
            done

            echo "ERROR: \$target failed download and/or validation after ${max_download_attempts} download attempts. See \$aria2_log" >&2
            rm -f "\$target" "\$target.aria2"
            record_failed_download "\$run_id" "\$fastq_file" "DOWNLOAD_OR_VALIDATION_FAILED_AFTER_RETRIES"
            return 1
        }

        # 下载 FASTQ；只有 gzip -t 和 ENA fastq_md5 都通过才写入下游清单
        if [[ "\$layout" == "PE" ]]; then
            for i in 0 1; do
                url="https://\${urls[\$i]}"
                target="fastqs/\${run_id}_\$((i+1)).fastq.gz"

                if ! download_and_verify_fastq "\$url" "\$target" "\${mds[\$i]}" "\$run_id"; then
                    echo "ERROR: Download/validation failed for \$target; failing this sample so Nextflow can retry or ignore it." >&2
                    exit 1
                fi
            done

            echo -e "${meta.id}\\t${species_name}\\t${ref_name}\\t\$run_id\\tPE\\t\$PWD/fastqs/\${run_id}_1.fastq.gz\\t\$PWD/fastqs/\${run_id}_2.fastq.gz\\tPLACEHOLDER_EXPECTED_RUNS" >> fastq_pairs.tmp
        elif [[ "\$layout" == "SE" ]]; then
            url="https://\${urls[0]}"
            target="fastqs/\${run_id}.fastq.gz"

            if ! download_and_verify_fastq "\$url" "\$target" "\${mds[0]}" "\$run_id"; then
                echo "ERROR: Download/validation failed for \$target; failing this sample so Nextflow can retry or ignore it." >&2
                exit 1
            fi

            echo -e "${meta.id}\\t${species_name}\\t${ref_name}\\t\$run_id\\tSE\\t\$PWD/fastqs/\${run_id}.fastq.gz\\t.\\tPLACEHOLDER_EXPECTED_RUNS" >> fastq_pairs.tmp
        fi
    done

    retained_count=\$(wc -l < fastq_pairs.tmp | awk '{print \$1}')

    if [[ "\$retained_count" -eq 0 ]]; then
        echo "ERROR: No usable FASTQ runs found for \${acc}." >&2
        exit 1
    fi

    echo -e "sample_id\\tspecies_name\\tref_name\\trun_id\\tlayout\\tr1\\tr2\\texpected_pairs" > fastq_pairs.tsv
    awk -v n="\$retained_count" 'BEGIN{FS=OFS="\t"} { \$8=n; print }' fastq_pairs.tmp >> fastq_pairs.tsv
    rm -f fastq_pairs.tmp

    trap - INT TERM HUP QUIT ERR
    """
}

process ALIGN_AND_SORT {
    tag "${meta.id}.${meta.pair_id}"
    label 'alignment_related'

    // 建议在 config 中将 exitStatus 1 加入 retry 策略
    // errorStrategy = { task.exitStatus in [1, 137, 140, 143] ? 'retry' : 'finish' }

    input:
    tuple val(meta), val(species_name), val(ref_name), path(reads)

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}.${meta.pair_id}.sorted.bam"), path("${meta.id}.${meta.pair_id}.sorted.bam.bai"), emit: bam

    script:
    def ref_file = "${params.global_ref_dir}/${ref_name}.fa"
    def bam_output = "${meta.id}.${meta.pair_id}.sorted.bam"
    def bwa_inputs = meta.layout == 'PE' ? "\"${reads[0]}\" \"${reads[1]}\"" : "\"${reads[0]}\""
    """
    #!/usr/bin/env bash
    set -Eeuo pipefail

    cleanup_alignment_outputs() {
        echo "WARN: ALIGN_AND_SORT interrupted or failed for ${meta.id}.${meta.pair_id}; removing partial BAM outputs." >&2
        rm -f "${bam_output}" "${bam_output}.bai" "${bam_output}.csi"
        rm -rf "tmp_sort_${meta.pair_id}"
    }

    trap cleanup_alignment_outputs INT TERM HUP QUIT ERR

    # Clean stale partial outputs from a previous failed/cancelled attempt before retrying.
    rm -f "${bam_output}" "${bam_output}.bai" "${bam_output}.csi"
    rm -rf "tmp_sort_${meta.pair_id}"
    
    echo "INFO: FASTQ layout: ${meta.layout}"
    echo "INFO: Validating input FASTQ integrity..."
    # 快速检查 gzip 文件是否损坏，如果是损坏文件则在此处提前终止报错
    if [[ "${meta.layout}" == "PE" ]]; then
        gzip -t "${reads[0]}" "${reads[1]}"
    elif [[ "${meta.layout}" == "SE" ]]; then
        gzip -t "${reads[0]}"
    else
        echo "ERROR: Unsupported FASTQ layout: ${meta.layout}" >&2
        exit 1
    fi

    for ext in amb ann bwt pac sa; do
        if [[ ! -s "${ref_file}.\${ext}" ]]; then
            echo "ERROR: Missing BWA index: ${ref_file}.\${ext}" >&2
            echo "Please run: bwa index ${ref_file}" >&2
            exit 1
        fi
    done

    if [[ ! -s "${ref_file}.fai" ]]; then
        echo "ERROR: Missing samtools faidx: ${ref_file}.fai" >&2
        echo "Please run: samtools faidx ${ref_file}" >&2
        exit 1
    fi

    # 创建独立的临时目录，避免多任务并发时可能产生的文件冲突
    # 使用当前工作目录下的 tmp，确保存储空间足够（通常比 /tmp 大）
    mkdir -p "tmp_sort_${meta.pair_id}"
    TMP_DIR="tmp_sort_${meta.pair_id}"

    echo "INFO: Starting BWA alignment and Samtools sort..."
    echo "INFO: Samtools sort threads: ${params.align_sort_threads}"
    echo "INFO: Per-thread sort memory: ${params.align_sort_mem_per_thread}"

    # ===== 执行核心管道 =====
    bwa mem -K 100000000 -v 3 -t ${task.cpus} -M -Y \\
      -R "@RG\\tID:${meta.id}.${meta.pair_id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}" \\
      "${ref_file}" ${bwa_inputs} | \\
    samtools sort -@ ${params.align_sort_threads} -m ${params.align_sort_mem_per_thread} \\
      -T "\${TMP_DIR}/sort_prefix" \\
      -o "${bam_output}" -

    echo "INFO: Indexing BAM..."
    samtools index -@ ${task.cpus} "${bam_output}"
    
    echo "INFO: Quick-checking BAM integrity..."
    samtools quickcheck "${bam_output}"

    # 成功完成后清理临时目录
    rm -rf "\${TMP_DIR}"

    trap - INT TERM HUP QUIT ERR
    """
}

process ALIGN_LARGE_FASTQ_STREAMING_CHUNKS {
    tag "${meta.id}.${meta.pair_id}"
    label 'alignment_related'

    input:
    tuple val(meta), val(species_name), val(ref_name), path(reads)

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}.${meta.pair_id}.sorted.bam"), path("${meta.id}.${meta.pair_id}.sorted.bam.bai"), emit: bam

    script:
    def ref_file = "${params.global_ref_dir}/${ref_name}.fa"
    def r1 = reads[0]
    def r2 = meta.layout == 'PE' ? reads[1] : null
    def cleanup_chunks = params.cleanup_chunk_intermediates ? 'true' : 'false'
    """
    #!/usr/bin/env bash
    set -Eeuo pipefail

    cleanup_streaming_outputs() {
        echo "WARN: ALIGN_LARGE_FASTQ_STREAMING_CHUNKS interrupted or failed for ${meta.id}.${meta.pair_id}; removing partial outputs." >&2
        rm -f "${meta.id}.${meta.pair_id}.sorted.bam" "${meta.id}.${meta.pair_id}.sorted.bam.bai" merged.tmp.bam merged.tmp.bam.bai chunk_bam.list
        rm -rf chunks chunk_bams sort_tmp
    }

    trap cleanup_streaming_outputs INT TERM HUP QUIT ERR

    rm -f "${meta.id}.${meta.pair_id}.sorted.bam" "${meta.id}.${meta.pair_id}.sorted.bam.bai" merged.tmp.bam merged.tmp.bam.bai chunk_bam.list
    rm -rf chunks chunk_bams sort_tmp
    mkdir -p chunks chunk_bams sort_tmp

    echo "INFO: ALIGN_LARGE_FASTQ_STREAMING_CHUNKS for ${meta.id}.${meta.pair_id}"
    echo "INFO: layout=${meta.layout}"
    echo "INFO: fastq_chunk_reads=${params.fastq_chunk_reads}"
    echo "INFO: align_sort_threads=${params.align_sort_threads}"
    echo "INFO: align_sort_mem_per_thread=${params.align_sort_mem_per_thread}"
    echo "INFO: cleanup_chunk_intermediates=${params.cleanup_chunk_intermediates}"
    echo "INFO: SLURM_JOB_ID=\${SLURM_JOB_ID:-}"
    echo "INFO: HOSTNAME=\$(hostname)"
    df -h .
    df -ih .
    ulimit -a
    ls -lh ${r1} ${r2 ?: ''}

    if command -v pigz >/dev/null 2>&1; then
        GZIP_CMD="pigz -p ${task.cpus}"
        ZCAT_CMD="pigz -dc"
    else
        GZIP_CMD="gzip"
        ZCAT_CMD="gzip -dc"
    fi
    export GZIP_CMD ZCAT_CMD
    echo "INFO: compression command: \${GZIP_CMD}"
    echo "INFO: decompression command: \${ZCAT_CMD}"

    if [[ "${meta.layout}" != "PE" && "${meta.layout}" != "SE" ]]; then
        echo "ERROR: Unsupported FASTQ layout: ${meta.layout}" >&2
        exit 1
    fi

    for ext in amb ann bwt pac sa; do
        [[ -s "${ref_file}.\${ext}" ]] || { echo "ERROR: Missing BWA index: ${ref_file}.\${ext}" >&2; exit 1; }
    done
    [[ -s "${ref_file}.fai" ]] || { echo "ERROR: Missing samtools faidx: ${ref_file}.fai" >&2; exit 1; }

    export SAMPLE_ID="${meta.id}"
    export PAIR_ID="${meta.pair_id}"
    export LAYOUT="${meta.layout}"
    export R1_FASTQ="${r1}"
    export R2_FASTQ="${r2 ?: ''}"
    export REF_FILE="${ref_file}"
    export CHUNK_READS="${params.fastq_chunk_reads}"
    export BWA_THREADS="${task.cpus}"
    export SORT_THREADS="${params.align_sort_threads}"
    export SORT_MEM="${params.align_sort_mem_per_thread}"
    export CLEANUP_CHUNK_INTERMEDIATES="${cleanup_chunks}"

    python <<'PY_STREAM_ALIGN'
import gzip
import os
import shlex
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

sample_id = os.environ["SAMPLE_ID"]
pair_id = os.environ["PAIR_ID"]
layout = os.environ["LAYOUT"]
r1_fastq = os.environ["R1_FASTQ"]
r2_fastq = os.environ.get("R2_FASTQ", "")
ref_file = os.environ["REF_FILE"]
chunk_reads = int(os.environ["CHUNK_READS"])
bwa_threads = os.environ["BWA_THREADS"]
sort_threads = os.environ["SORT_THREADS"]
sort_mem = os.environ["SORT_MEM"]
gzip_cmd = shlex.split(os.environ["GZIP_CMD"])
cleanup = os.environ.get("CLEANUP_CHUNK_INTERMEDIATES", "true").lower() == "true"

Path("chunks").mkdir(exist_ok=True)
Path("chunk_bams").mkdir(exist_ok=True)
Path("sort_tmp").mkdir(exist_ok=True)


def log(msg):
    print("{} INFO: {}".format(datetime.now().isoformat(timespec="seconds"), msg), flush=True)


def read_record(handle, label):
    rec = [handle.readline() for _ in range(4)]
    if rec[0] == "":
        if any(line != "" for line in rec[1:]):
            raise RuntimeError("Truncated FASTQ record in {}".format(label))
        return None
    if any(line == "" for line in rec):
        raise RuntimeError("Truncated FASTQ record in {}".format(label))
    return rec


def compress_fastq(path):
    subprocess.run(gzip_cmd + [str(path)], check=True)
    return Path(str(path) + ".gz")


def run_cmd(cmd):
    log("Running: " + cmd)
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)


def align_chunk(chunk_id, r1_gz, r2_gz=None):
    bam = Path("chunk_bams") / "{}.{}.{}.sorted.bam".format(sample_id, pair_id, chunk_id)
    tmp_prefix = Path("sort_tmp") / "{}.{}".format(pair_id, chunk_id)
    bs = chr(92)
    escaped_tab = bs + "t"
    literal_tab = chr(9)
    rg = (
        "@RG"
        + escaped_tab + "ID:" + sample_id + "." + pair_id
        + escaped_tab + "SM:" + sample_id
        + escaped_tab + "PL:ILLUMINA"
        + escaped_tab + "LB:" + sample_id
    )
    if literal_tab in rg:
        raise RuntimeError("RG string contains literal tab characters; expected escaped backslash-t sequences")
    if escaped_tab not in rg:
        raise RuntimeError("RG string does not contain escaped backslash-t sequences")
    log("BWA RG string repr: " + repr(rg))
    log("BWA RG string: " + rg)
    inputs = [shlex.quote(str(r1_gz))]
    if r2_gz is not None:
        inputs.append(shlex.quote(str(r2_gz)))
    cmd = (
        "bwa mem -K 100000000 -v 3 -t {threads} -M -Y -R {rg} {ref} {inputs} | "
        "samtools sort -@ {sort_threads} -m {sort_mem} -T {tmp} -o {bam} -"
    ).format(
        threads=shlex.quote(str(bwa_threads)),
        rg=shlex.quote(rg),
        ref=shlex.quote(ref_file),
        inputs=" ".join(inputs),
        sort_threads=shlex.quote(str(sort_threads)),
        sort_mem=shlex.quote(str(sort_mem)),
        tmp=shlex.quote(str(tmp_prefix)),
        bam=shlex.quote(str(bam)),
    )
    run_cmd(cmd)
    run_cmd("samtools index {}".format(shlex.quote(str(bam))))
    run_cmd("samtools quickcheck -v {}".format(shlex.quote(str(bam))))
    log("chunk {} BAM size: {} bytes".format(chunk_id, bam.stat().st_size))
    if cleanup:
        for fq in [r1_gz, r2_gz]:
            if fq:
                Path(fq).unlink(missing_ok=True)
        for candidate in Path("sort_tmp").glob("{}.{}*".format(pair_id, chunk_id)):
            if candidate.is_dir():
                shutil.rmtree(candidate)
            else:
                candidate.unlink(missing_ok=True)
    subprocess.run(["du", "-sh", "."], check=False)
    return bam


def write_and_align_pe():
    bams = []
    chunk_index = 1
    with gzip.open(r1_fastq, "rt") as r1, gzip.open(r2_fastq, "rt") as r2:
        while True:
            chunk_id = "chunk_{:06d}".format(chunk_index)
            r1_path = Path("chunks") / "{}_R1.fastq".format(chunk_id)
            r2_path = Path("chunks") / "{}_R2.fastq".format(chunk_id)
            pairs = 0
            log("chunk {} start".format(chunk_id))
            with r1_path.open("wt") as r1_out, r2_path.open("wt") as r2_out:
                for _ in range(chunk_reads):
                    rec1 = read_record(r1, "R1")
                    rec2 = read_record(r2, "R2")
                    if rec1 is None and rec2 is None:
                        break
                    if rec1 is None or rec2 is None:
                        raise RuntimeError("R1 and R2 have different record counts")
                    r1_out.writelines(rec1)
                    r2_out.writelines(rec2)
                    pairs += 1
            if pairs == 0:
                r1_path.unlink(missing_ok=True)
                r2_path.unlink(missing_ok=True)
                break
            log("chunk {} read pairs written: {}".format(chunk_id, pairs))
            r1_gz = compress_fastq(r1_path)
            r2_gz = compress_fastq(r2_path)
            log("chunk {} FASTQ sizes: R1={} bytes R2={} bytes".format(chunk_id, r1_gz.stat().st_size, r2_gz.stat().st_size))
            bams.append(align_chunk(chunk_id, r1_gz, r2_gz))
            log("chunk {} end".format(chunk_id))
            chunk_index += 1
    return bams


def write_and_align_se():
    bams = []
    chunk_index = 1
    with gzip.open(r1_fastq, "rt") as r1:
        while True:
            chunk_id = "chunk_{:06d}".format(chunk_index)
            chunk_path = Path("chunks") / "{}.fastq".format(chunk_id)
            reads = 0
            log("chunk {} start".format(chunk_id))
            with chunk_path.open("wt") as out:
                for _ in range(chunk_reads):
                    rec = read_record(r1, "SE")
                    if rec is None:
                        break
                    out.writelines(rec)
                    reads += 1
            if reads == 0:
                chunk_path.unlink(missing_ok=True)
                break
            log("chunk {} reads written: {}".format(chunk_id, reads))
            chunk_gz = compress_fastq(chunk_path)
            log("chunk {} FASTQ size: {} bytes".format(chunk_id, chunk_gz.stat().st_size))
            bams.append(align_chunk(chunk_id, chunk_gz))
            log("chunk {} end".format(chunk_id))
            chunk_index += 1
    return bams

if layout == "PE":
    chunk_bams = write_and_align_pe()
elif layout == "SE":
    chunk_bams = write_and_align_se()
else:
    raise RuntimeError("Unsupported FASTQ layout: {}".format(layout))

if not chunk_bams:
    raise RuntimeError("No chunk BAMs produced")

with open("chunk_bam.list", "wt") as out:
    for bam in sorted(chunk_bams):
        print(str(bam), file=out)

log("Merging {} chunk BAMs".format(len(chunk_bams)))
run_cmd("samtools merge -@ {} -f -b chunk_bam.list merged.tmp.bam".format(shlex.quote(str(bwa_threads))))
run_cmd("samtools index merged.tmp.bam")
run_cmd("samtools quickcheck -v merged.tmp.bam")
Path("merged.tmp.bam").replace("{}.{}.sorted.bam".format(sample_id, pair_id))
Path("merged.tmp.bam.bai").replace("{}.{}.sorted.bam.bai".format(sample_id, pair_id))
run_cmd("samtools quickcheck -v {}".format(shlex.quote("{}.{}.sorted.bam".format(sample_id, pair_id))))

if cleanup:
    shutil.rmtree("chunks", ignore_errors=True)
    shutil.rmtree("sort_tmp", ignore_errors=True)
    for path in Path("chunk_bams").glob("*.sorted.bam*"):
        path.unlink(missing_ok=True)
    Path("chunk_bam.list").unlink(missing_ok=True)
    try:
        Path("chunk_bams").rmdir()
    except OSError:
        pass
PY_STREAM_ALIGN

    trap - INT TERM HUP QUIT ERR
    """
}

process SPLIT_FASTQ_CHUNKS {
    tag "${meta.id}.${meta.pair_id}"
    label 'alignment_related'

    input:
    tuple val(meta), val(species_name), val(ref_name), path(reads)

    output:
    path "chunks.tsv", emit: chunks_tsv
    path "chunks/*.fastq.gz", emit: chunk_fastqs

    script:
    def chunk_lines = (params.fastq_chunk_reads as long) * 4L
    def r1 = reads[0]
    def r2 = meta.layout == 'PE' ? reads[1] : null
    """
    #!/usr/bin/env bash
    set -Eeuo pipefail

    cleanup_split_outputs() {
        echo "WARN: SPLIT_FASTQ_CHUNKS interrupted or failed for ${meta.id}.${meta.pair_id}; removing partial chunks." >&2
        rm -rf chunks chunks.tsv chunk_ids.txt
    }

    trap cleanup_split_outputs INT TERM HUP QUIT ERR

    rm -rf chunks chunks.tsv chunk_ids.txt
    mkdir -p chunks

    echo "INFO: SPLIT_FASTQ_CHUNKS diagnostics for ${meta.id}.${meta.pair_id}"
    echo "INFO: SLURM_JOB_ID=\${SLURM_JOB_ID:-}"
    echo "INFO: HOSTNAME=\$(hostname)"
    df -h .
    df -ih .
    ulimit -a
    echo "INFO: Chunked alignment enabled: ${params.enable_chunked_alignment}"
    echo "INFO: FASTQ size threshold GB: ${params.chunked_alignment_size_threshold_gb}"
    echo "INFO: Total FASTQ size GB: ${meta.total_fastq_size_gb}"
    echo "INFO: Requested reads per chunk: ${params.fastq_chunk_reads}"
    echo "INFO: Lines per chunk: ${chunk_lines}"
    echo "INFO: FASTQ R1 size bytes: \$(stat -c%s '${r1}')"
    if [[ "${meta.layout}" == "PE" ]]; then
        echo "INFO: FASTQ R2 size bytes: \$(stat -c%s '${r2}')"
    fi

    if command -v pigz >/dev/null 2>&1; then
        GZIP_CMD="pigz -p ${task.cpus}"
    else
        GZIP_CMD="gzip"
    fi
    export GZIP_CMD

    echo "INFO: compression command: \${GZIP_CMD}"

    if [[ "${meta.layout}" == "PE" ]]; then
        zcat "${r1}" | split -d -a 6 --numeric-suffixes=1 -l ${chunk_lines} --filter='\${GZIP_CMD} > chunks/\${FILE}_R1.fastq.gz' - chunk_
        zcat "${r2}" | split -d -a 6 --numeric-suffixes=1 -l ${chunk_lines} --filter='\${GZIP_CMD} > chunks/\${FILE}_R2.fastq.gz' - chunk_

        for r1_chunk in chunks/chunk_*_R1.fastq.gz; do
            [[ -e "\${r1_chunk}" ]] || { echo "ERROR: No R1 chunks produced." >&2; exit 1; }
            chunk_id=\$(basename "\${r1_chunk}" _R1.fastq.gz)
            r2_chunk="chunks/\${chunk_id}_R2.fastq.gz"
            if [[ ! -s "\${r2_chunk}" ]]; then
                echo "ERROR: Missing paired R2 chunk for \${chunk_id}" >&2
                exit 1
            fi
            gzip -t "\${r1_chunk}" "\${r2_chunk}"
            echo "\${chunk_id}" >> chunk_ids.txt
        done
    elif [[ "${meta.layout}" == "SE" ]]; then
        zcat "${r1}" | split -d -a 6 --numeric-suffixes=1 -l ${chunk_lines} --filter='\${GZIP_CMD} > chunks/\${FILE}.fastq.gz' - chunk_

        for chunk in chunks/chunk_*.fastq.gz; do
            [[ -e "\${chunk}" ]] || { echo "ERROR: No SE chunks produced." >&2; exit 1; }
            chunk_id=\$(basename "\${chunk}" .fastq.gz)
            gzip -t "\${chunk}"
            echo "\${chunk_id}" >> chunk_ids.txt
        done
    else
        echo "ERROR: Unsupported FASTQ layout: ${meta.layout}" >&2
        exit 1
    fi

    sort -u chunk_ids.txt -o chunk_ids.txt
    expected_chunks=\$(wc -l < chunk_ids.txt | awk '{print \$1}')
    echo "INFO: Number of chunks: \${expected_chunks}"
    echo "INFO: Chunk IDs: \$(paste -sd, chunk_ids.txt)"

    echo -e "sample_id\\tspecies_name\\tref_name\\trun_id\\tlayout\\tchunk_id\\tr1\\tr2\\texpected_chunks\\texpected_pairs" > chunks.tsv
    while read -r chunk_id; do
        if [[ "${meta.layout}" == "PE" ]]; then
            echo -e "${meta.id}\\t${species_name}\\t${ref_name}\\t${meta.pair_id}\\t${meta.layout}\\t\${chunk_id}\\t\$PWD/chunks/\${chunk_id}_R1.fastq.gz\\t\$PWD/chunks/\${chunk_id}_R2.fastq.gz\\t\${expected_chunks}\\t${meta.n_pairs}" >> chunks.tsv
        else
            echo -e "${meta.id}\\t${species_name}\\t${ref_name}\\t${meta.pair_id}\\t${meta.layout}\\t\${chunk_id}\\t\$PWD/chunks/\${chunk_id}.fastq.gz\\t.\\t\${expected_chunks}\\t${meta.n_pairs}" >> chunks.tsv
        fi
    done < chunk_ids.txt

    trap - INT TERM HUP QUIT ERR
    """
}

process ALIGN_AND_SORT_CHUNK {
    tag "${meta.id}.${meta.pair_id}.${meta.chunk_id}"
    label 'alignment_related'

    input:
    tuple val(meta), val(species_name), val(ref_name), path(reads)

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}.${meta.pair_id}.${meta.chunk_id}.sorted.bam"), path("${meta.id}.${meta.pair_id}.${meta.chunk_id}.sorted.bam.bai"), emit: bam

    script:
    def ref_file = "${params.global_ref_dir}/${ref_name}.fa"
    def bam_output = "${meta.id}.${meta.pair_id}.${meta.chunk_id}.sorted.bam"
    def bwa_inputs = meta.layout == 'PE' ? "\"${reads[0]}\" \"${reads[1]}\"" : "\"${reads[0]}\""
    """
    #!/usr/bin/env bash
    set -Eeuo pipefail

    cleanup_chunk_alignment_outputs() {
        echo "WARN: ALIGN_AND_SORT_CHUNK interrupted or failed for ${meta.id}.${meta.pair_id}.${meta.chunk_id}; removing partial BAM outputs." >&2
        rm -f "${bam_output}" "${bam_output}.bai" "${bam_output}.csi"
        rm -rf "tmp_sort_${meta.id}_${meta.pair_id}_${meta.chunk_id}"
    }

    trap cleanup_chunk_alignment_outputs INT TERM HUP QUIT ERR

    rm -f "${bam_output}" "${bam_output}.bai" "${bam_output}.csi"
    rm -rf "tmp_sort_${meta.id}_${meta.pair_id}_${meta.chunk_id}"

    echo "INFO: ALIGN_AND_SORT_CHUNK diagnostics for ${meta.id}.${meta.pair_id}.${meta.chunk_id}"
    echo "INFO: SLURM_JOB_ID=\${SLURM_JOB_ID:-}"
    echo "INFO: HOSTNAME=\$(hostname)"
    df -h .
    df -ih .
    ulimit -a
    echo "INFO: FASTQ layout: ${meta.layout}"
    echo "INFO: Chunk ID: ${meta.chunk_id}"
    echo "INFO: Expected chunks for run: ${meta.n_chunks}"

    if [[ "${meta.layout}" == "PE" ]]; then
        gzip -t "${reads[0]}" "${reads[1]}"
    else
        gzip -t "${reads[0]}"
    fi

    for ext in amb ann bwt pac sa; do
        [[ -s "${ref_file}.\${ext}" ]] || { echo "ERROR: Missing BWA index: ${ref_file}.\${ext}" >&2; exit 1; }
    done
    [[ -s "${ref_file}.fai" ]] || { echo "ERROR: Missing samtools faidx: ${ref_file}.fai" >&2; exit 1; }

    mkdir -p "tmp_sort_${meta.id}_${meta.pair_id}_${meta.chunk_id}"

    bwa mem -K 100000000 -v 3 -t ${task.cpus} -M -Y \\
      -R "@RG\\tID:${meta.id}.${meta.pair_id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}" \\
      "${ref_file}" ${bwa_inputs} | \\
    samtools sort -@ ${params.align_sort_threads} -m ${params.align_sort_mem_per_thread} \\
      -T "tmp_sort_${meta.id}_${meta.pair_id}_${meta.chunk_id}/sort_prefix" \\
      -o "${bam_output}" -

    samtools index -@ ${task.cpus} "${bam_output}"
    samtools quickcheck "${bam_output}"
    rm -rf "tmp_sort_${meta.id}_${meta.pair_id}_${meta.chunk_id}"

    trap - INT TERM HUP QUIT ERR
    """
}

process MERGE_RUN_CHUNK_BAMS {
    tag "${meta.id}.${meta.pair_id}"
    label 'merge_related'

    input:
    tuple val(meta), val(species_name), val(ref_name), path(bams)

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}.${meta.pair_id}.sorted.bam"), path("${meta.id}.${meta.pair_id}.sorted.bam.bai"), emit: bam

    script:
    def bam_list = bams.join(' ')
    def merged_name = "${meta.id}.${meta.pair_id}.sorted.bam"
    def tmp_name = "${meta.id}.${meta.pair_id}.merged.tmp.bam"
    """
    #!/usr/bin/env bash
    set -Eeuo pipefail

    cleanup_run_merge_outputs() {
        echo "WARN: MERGE_RUN_CHUNK_BAMS interrupted or failed for ${meta.id}.${meta.pair_id}; removing partial outputs." >&2
        rm -f "${tmp_name}" "${tmp_name}.bai" "${tmp_name}.csi"
        rm -f "${merged_name}" "${merged_name}.bai" "${merged_name}.csi"
    }

    trap cleanup_run_merge_outputs INT TERM HUP QUIT ERR

    rm -f "${tmp_name}" "${tmp_name}.bai" "${tmp_name}.csi"
    rm -f "${merged_name}" "${merged_name}.bai" "${merged_name}.csi"

    echo "INFO: MERGE_RUN_CHUNK_BAMS diagnostics for ${meta.id}.${meta.pair_id}"
    echo "INFO: SLURM_JOB_ID=\${SLURM_JOB_ID:-}"
    echo "INFO: HOSTNAME=\$(hostname)"
    df -h .
    df -ih .
    ulimit -a
    echo "INFO: Number of chunk BAMs: ${bams.size()}"
    echo "INFO: Chunk BAMs: ${bam_list}"

    if [[ ${bams.size()} -ne ${meta.n_chunks} ]]; then
        echo "ERROR: expected ${meta.n_chunks} chunk BAMs but received ${bams.size()}" >&2
        exit 1
    fi
    printf '%s\n' ${bam_list} > chunk_bams.list
    echo "INFO: Ordered chunk BAM list:"
    cat chunk_bams.list

    samtools merge -@ ${task.cpus} -f -b chunk_bams.list "${tmp_name}"
    samtools index -@ ${task.cpus} "${tmp_name}"
    samtools quickcheck "${tmp_name}"

    mv "${tmp_name}" "${merged_name}"
    mv "${tmp_name}.bai" "${merged_name}.bai"
    samtools quickcheck "${merged_name}"

    trap - INT TERM HUP QUIT ERR
    """
}

process MERGE_BAMS {
    tag "${meta.id}"
    label 'merge_related'

    input:
    tuple val(meta), val(species_name), val(ref_name), path(bams)

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}.merged.bam"), path("${meta.id}.merged.bam.bai"), emit: merged_bam

    script:
    def bam_list = bams.join(' ')
    def merged_name = "${meta.id}.merged.bam"
    def tmp_name = "${meta.id}.merged.tmp.bam"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    cleanup_partial_outputs() {
        echo "WARN: MERGE_BAMS interrupted or failed for ${meta.id}; removing partial outputs." >&2
        rm -f "${tmp_name}" "${tmp_name}.bai" "${tmp_name}.csi"
        rm -f "${merged_name}" "${merged_name}.bai" "${merged_name}.csi"
    }

    trap cleanup_partial_outputs INT TERM HUP QUIT ERR

    # Clean stale partial outputs from previous cancelled attempts.
    rm -f "${tmp_name}" "${tmp_name}.bai" "${tmp_name}.csi"
    rm -f "${merged_name}" "${merged_name}.bai" "${merged_name}.csi"

    echo "INFO: Merging BAMs for ${meta.id} into temporary file ${tmp_name}"
    samtools merge -f -@ ${task.cpus} "${tmp_name}" ${bam_list}

    echo "INFO: Indexing temporary merged BAM"
    samtools index -@ ${task.cpus} "${tmp_name}"

    echo "INFO: Quick-checking temporary merged BAM"
    samtools quickcheck "${tmp_name}"

    echo "INFO: Promoting temporary merged BAM to final output"
    mv "${tmp_name}" "${merged_name}"
    mv "${tmp_name}.bai" "${merged_name}.bai"

    echo "INFO: Final quickcheck"
    samtools quickcheck "${merged_name}"

    trap - INT TERM HUP QUIT ERR
    """
}

process BAM_TO_CRAM {
    tag "${meta.id}"
    label 'alignment_related'
    // The marker is created last and is the publication-completion signal. The
    // startup validator still waits for stable destination sizes because copy
    // publication is not an atomic rename on every shared filesystem.
    publishDir "${params.outdir}/${meta.id}/alignment", mode: 'copy', pattern: "*.{cram,crai,complete}"

    input:
    tuple val(meta), val(species_name), val(ref_name), path(bam), path(bai)

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}.cram"), path("${meta.id}.cram.crai"), path("${meta.id}.cram.complete"), emit: cram

    script:
    def ref_file = "${params.global_ref_dir}/${ref_name}.fa"
    def cram_out = "${meta.id}.cram"
    def crai_out = "${meta.id}.cram.crai"
    
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # 1. 检查参考序列
    if [[ ! -f "${ref_file}" ]]; then
        echo "ERROR: Reference ${ref_file} not found" >&2
        exit 1
    fi

    # 2. BAM 转 CRAM
    # -C 选项代表输出 CRAM
    samtools view -@ ${task.cpus} -T "${ref_file}" -C \\
      -o "${cram_out}" "${bam}"
    
    # 3. 显式指定索引文件名输出
    # 这样可以确保输出文件名绝对符合 output 定义的 ${meta.id}.cram.crai
    samtools index -@ ${task.cpus} "${cram_out}" "${crai_out}"
    
    # 4. 验证 CRAM 文件是否损坏
    samtools quickcheck -v "${cram_out}"
    samtools idxstats --reference "${ref_file}" "${cram_out}" >/dev/null

    cram_size=\$(stat -c%s "${cram_out}")
    crai_size=\$(stat -c%s "${crai_out}")
    samtools_version=\$(samtools --version | head -n1)
    cat > "${meta.id}.cram.complete.tmp" <<EOF
sample_id=${meta.id}
ref_name=${ref_name}
cram=${cram_out}
crai=${crai_out}
cram_size=\${cram_size}
crai_size=\${crai_size}
samtools_version=\${samtools_version}
completed_at=\$(date --iso-8601=seconds)
EOF
    mv "${meta.id}.cram.complete.tmp" "${meta.id}.cram.complete"
    """
}

workflow.onComplete {
    if (workflow.success) {
        log.info "Pipeline completed successfully at: ${workflow.complete}"
    }
}

workflow.onError {
    log.error "Pipeline failed! Check the log file for errors."
}
