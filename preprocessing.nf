#!/usr/bin/env nextflow
nextflow.enable.dsl=2

if (!params.containsKey('skip_existing_cram') || params.skip_existing_cram == null) {
    params.skip_existing_cram = true
}

def existingCramIsComplete = { String sampleId ->
    def cram = new File("${params.outdir}/${sampleId}/alignment/${sampleId}.cram")
    def crai = new File("${params.outdir}/${sampleId}/alignment/${sampleId}.cram.crai")

    if (!cram.isFile() || cram.length() == 0 || !crai.isFile() || crai.length() == 0) {
        return false
    }

    try {
        def proc = new ProcessBuilder('samtools', 'quickcheck', cram.absolutePath)
            .redirectErrorStream(true)
            .start()
        proc.waitFor()
        return proc.exitValue() == 0
    } catch (IOException e) {
        log.warn "Could not run 'samtools quickcheck' for existing CRAM ${cram}; will reprocess sample ${sampleId}."
        return false
    }
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
            if (params.skip_existing_cram && existingCramIsComplete(meta.id)) {
                log.info "Skipping sample ${meta.id}: complete CRAM and CRAI already exist under ${params.outdir}/${meta.id}/alignment"
                return false
            }
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
            def r1 = file(row.r1)
            def r2 = row.layout == 'PE' ? file(row.r2) : file(row.r1)
            tuple(meta, row.species_name, row.ref_name, r1, r2)
        }

    // 第三步：比对并排序
    ALIGN_AND_SORT(ch_fastq_pairs)

    // 第四步：按 Sample ID 分组，并强制校验数量
    ch_bam_grouped = ALIGN_AND_SORT.out.bam
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
    tuple val(meta), val(species_name), val(ref_name), path(r1), path(r2)

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}.${meta.pair_id}.sorted.bam"), path("${meta.id}.${meta.pair_id}.sorted.bam.bai"), emit: bam

    script:
    def ref_file = "${params.global_ref_dir}/${ref_name}.fa"
    // 将比例从 0.7 降至 0.6，预留更多 buffer 防止系统层面杀进程
    def sort_mem = task.memory ? "${(task.memory.toGiga() * 0.6 / task.cpus).toInteger()}G" : "2G"
    def bam_output = "${meta.id}.${meta.pair_id}.sorted.bam"
    def bwa_inputs = meta.layout == 'PE' ? "\"${r1}\" \"${r2}\"" : "\"${r1}\""
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
        gzip -t "${r1}" "${r2}"
    elif [[ "${meta.layout}" == "SE" ]]; then
        gzip -t "${r1}"
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
    echo "INFO: Per-thread sort memory: ${sort_mem}"

    # ===== 执行核心管道 =====
    bwa mem -K 100000000 -v 3 -t ${task.cpus} -M -Y \\
      -R "@RG\\tID:${meta.id}.${meta.pair_id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}" \\
      "${ref_file}" ${bwa_inputs} | \\
    samtools sort -@ ${task.cpus} -m ${sort_mem} \\
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
    publishDir "${params.outdir}/${meta.id}/alignment", mode: 'copy', pattern: "*.{cram,crai}"

    input:
    tuple val(meta), val(species_name), val(ref_name), path(bam), path(bai)

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}.cram"), path("${meta.id}.cram.crai"), emit: cram

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
    samtools quickcheck "${cram_out}"
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
