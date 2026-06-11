#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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
======================================
"""

// 1. 从 TSV 读取初始样本信息
ch_samples = Channel.fromPath(params.sample_tsv)
    .splitCsv(header: false, sep: '\t')
    .filter { row -> row.size() >= 3 && row[0]?.trim() }
    .map { row ->
        def meta = [id: row[0].trim()]
        def species = row[1].trim()
        def ref_name = row[2].trim()
        tuple(meta, species, ref_name)
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
                n_pairs : row.expected_pairs.toInteger()
            ]
            // 必须使用 file() 包装路径，Nextflow 才能在进程间正确传递文件
            def r1 = file(row.r1)
            def r2 = file(row.r2)
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
    
    // 建议增加重试策略，应对 EBI 服务器偶尔的拒绝访问
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(meta), val(species_name), val(ref_name)

    output:
    path "fastq_pairs.tsv", emit: pairs_tsv
    path "fastqs/*.fastq.gz", emit: fastq_files

    script:
    // 将你的 aria2c 路径定义为 Nextflow 变量
    def aria2_bin = "/home/lt692/.conda/envs/aria2_env/bin/aria2c"
    
    """
    #!/usr/bin/env bash
    set -euo pipefail

    mkdir -p fastqs
    ena_base="https://www.ebi.ac.uk/ena/portal/api"
    acc="${meta.id}"

    # 1. 获取 Run IDs
    url_report="\${ena_base}/filereport?accession=\${acc}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5&format=tsv"
    curl -fsSL "\$url_report" -o report.tsv
    
    if [[ \$(tail -n +2 report.tsv | wc -l) -eq 0 ]]; then
        echo "ERROR: No runs found for \${acc}" >&2; exit 1
    fi

    # 2. 准备输出清单
    echo -e "sample_id\\tspecies_name\\tref_name\\trun_id\\tr1\\tr2\\texpected_pairs" > fastq_pairs.tsv
    valid_count=\$(tail -n +2 report.tsv | wc -l)

    # 3. 循环下载并验证 MD5
    tail -n +2 report.tsv | while IFS=\$'\\t' read -r run_id ftp_urls md5s; do
        IFS=';' read -r -a urls <<< "\$ftp_urls"
        IFS=';' read -r -a mds <<< "\$md5s"

        if [[ \${#urls[@]} -ne 2 ]]; then
            echo "ERROR: Run \$run_id is not paired-end." >&2; exit 1
        fi

        # 下载 R1 和 R2
        for i in 0 1; do
            url="https://\${urls[\$i]}"
            target="fastqs/\${run_id}_\$((i+1)).fastq.gz"
            
            echo "Downloading \$url with aria2c..."
            
            # 使用 aria2c 增强下载
            # -x 10 -s 10: 10个线程加速
            # -c: 核心！支持断点续传
            # --summary-interval=0: 减少日志刷屏
            ${aria2_bin} -x 10 -s 10 -c -m 0 --retry-wait 5 \\
                --summary-interval=0 -d fastqs -o "\${run_id}_\$((i+1)).fastq.gz" "\$url"

            # 验证 MD5，如果失败则删除文件并报错退出
            if ! echo "\${mds[\$i]}  \$target" | md5sum -c -; then
                echo "ERROR: MD5 check failed for \$target. Removing corrupted file." >&2
                rm -f "\$target"
                exit 1
            fi
        done

        # 写入清单
        echo -e "${meta.id}\\t${species_name}\\t${ref_name}\\t\$run_id\\t\$PWD/fastqs/\${run_id}_1.fastq.gz\\t\$PWD/fastqs/\${run_id}_2.fastq.gz\\t\$valid_count" >> fastq_pairs.tsv
    done
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
    tuple val(meta), val(species_name), val(ref_name), path("*.bam"), path("*.bai"), emit: bam

    script:
    def ref_file = "${params.global_ref_dir}/${ref_name}.fasta"
    // 将比例从 0.7 降至 0.6，预留更多 buffer 防止系统层面杀进程
    def sort_mem = task.memory ? "${(task.memory.toGiga() * 0.6 / task.cpus).toInteger()}G" : "2G"
    def bam_output = "${meta.id}.${meta.pair_id}.sorted.bam"
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    echo "INFO: Validating input FASTQ integrity..."
    # 快速检查 gzip 文件是否损坏，如果是损坏文件则在此处提前终止报错
    gzip -t "${r1}" "${r2}"

    # 创建独立的临时目录，避免多任务并发时可能产生的文件冲突
    # 使用当前工作目录下的 tmp，确保存储空间足够（通常比 /tmp 大）
    mkdir -p "tmp_sort_${meta.pair_id}"
    TMP_DIR="tmp_sort_${meta.pair_id}"

    echo "INFO: Starting BWA alignment and Samtools sort..."
    echo "INFO: Per-thread sort memory: ${sort_mem}"

    # ===== 执行核心管道 =====
    bwa mem -K 100000000 -v 3 -t ${task.cpus} -M -Y \\
      -R "@RG\\tID:${meta.id}.${meta.pair_id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}" \\
      "${ref_file}" "${r1}" "${r2}" | \\
    samtools sort -@ ${task.cpus} -m ${sort_mem} \\
      -T "\${TMP_DIR}/sort_prefix" \\
      -o "${bam_output}" -

    echo "INFO: Indexing BAM..."
    samtools index -@ ${task.cpus} "${bam_output}"
    
    echo "INFO: Quick-checking BAM integrity..."
    samtools quickcheck "${bam_output}"

    # 成功完成后清理临时目录
    rm -rf "\${TMP_DIR}"
    """
}

process MERGE_BAMS {
    tag "${meta.id}"
    label 'merge_related'

    input:
    tuple val(meta), val(species_name), val(ref_name), path(bams)

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}.merged.bam"), path("*.bai"), emit: merged_bam

    script:
    def bam_list = bams.join(' ')
    def merged_name = "${meta.id}.merged.bam"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    samtools merge -f -@ ${task.cpus} "${merged_name}" ${bam_list}
    samtools index -@ ${task.cpus} "${merged_name}"
    
    samtools quickcheck "${merged_name}"
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
    def ref_file = "${params.global_ref_dir}/${ref_name}.fasta"
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
    log.info "Pipeline completed successfully at: \${workflow.complete}"
}

workflow.onError {
    log.error "Pipeline failed! Check the log file for errors."
}