#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ====================================================================================
//                                  PIPELINE HEADER
// ====================================================================================
log.info """
PRIMATE MITOCHONDRIAL PIPELINE START
=================================
Sample TSV:          ${params.sample_tsv}
Output Directory:    ${params.outdir}
Mito Reference Dir:  ${params.ref_dir}
WDL Script:          ${params.wdl_script}
=================================
"""

// ====================================================================================
//                                  INPUT CHANNELS
// ====================================================================================
ch_samples = Channel.fromPath(params.sample_tsv)
    .splitCsv(header: false, sep: '\t')
    .filter { row -> row.size() >= 3 && row[0]?.trim() && row[1]?.trim() && row[2]?.trim() }
    .map { row ->
        def meta = [id: row[0].trim()]
        def species = row[1].trim()
        def ref_name = row[2].trim() // 第三列：参考序列名称
        tuple(meta, species, ref_name)
    }

ch_cromwell_conf = file("${baseDir}/cromwell.conf")

// ====================================================================================
//                                    WORKFLOW
// ====================================================================================
workflow {

    // --- Step 1: 若最终结果已存在则跳过 ---
    ch_samples
        .filter { meta, species, ref_name ->
            def cromwellOutputDir = file("${params.outdir}/${meta.id}/variant_calling/cromwell-executions")
            if (cromwellOutputDir.exists() && cromwellOutputDir.isDirectory()) {
                log.info "[SKIP] Sample '${meta.id}' already complete."
                return false
            }
            true
        }
        .set { ch_samples_to_process }

    // --- Step 2: 检查是否已有 CRAM，分流 ---
    def branched = ch_samples_to_process.branch { meta, species, ref_name ->
    def cram_file = file("${params.outdir}/${meta.id}/alignment/${meta.id}.cram")
    partial_run: cram_file.exists() && cram_file.isFile()
    full_run  : !(cram_file.exists() && cram_file.isFile())
  }

    // -------------------------------
// A) FULL RUN: download -> pairs
// -------------------------------

DOWNLOAD_FASTQ( branched.full_run )

DOWNLOAD_FASTQ.out
  .map { meta, species, ref_name, fastq_dir, done -> tuple(meta, species, ref_name, fastq_dir) }
  .flatMap { meta, species, ref_name, fastq_dir ->

      def r1s = file("${fastq_dir}/*_1.fastq*").collect().sort { it.name }
      def r2s = file("${fastq_dir}/*_2.fastq*").collect().sort { it.name }

      if (r1s.size() != r2s.size()) {
          log.warn "Sample ${meta.id}: number of R1 and R2 files mismatch: ${r1s.size()} vs ${r2s.size()}"
      }

      // 先匹配出“有效 pairs”
      def pairs = []
      r1s.each { r1 ->
          def stem = r1.name.replaceAll(/_1\.fastq.*/, '')
          def r2 = r2s.find { it.name.startsWith(stem) }
          if (r2) {
              pairs << [stem, r1, r2]
          } else {
              log.warn "Sample ${meta.id}: could not find R2 for ${r1.name}, skipping this pair"
          }
      }

      // 关键：给每个 pair_meta 写入 n_pairs（该 ERS 期望的 pair 数）
      def n_pairs = pairs.size()
      if (n_pairs == 0) {
          log.warn "[SKIP ERS] ${meta.id}: no valid FASTQ pairs found; skip."
          return []
      }

      pairs.collect { stem, r1, r2 ->
          def pair_meta = meta + [pair_id: stem, n_pairs: n_pairs]
          tuple(pair_meta, species, ref_name, r1, r2)
      }
  }
  .set { ch_fastq_pairs }

ALIGN_AND_SORT( ch_fastq_pairs )

// -----------------------------------------
// B) Collect BAMs per ERS and gate by n_pairs
//    - complete: merge -> cram
//    - incomplete: skip
// -----------------------------------------

def ch_bam_grouped = ALIGN_AND_SORT.out.bam
  .map { meta, species, ref_name, bam ->

      // meta is per pair; make sample-level meta2
      def meta2 = (meta instanceof Map) ? new LinkedHashMap(meta) : meta
      if (meta2 instanceof Map) meta2.remove('pair_id')

      // robust expected number
      def expected = (meta instanceof Map ? meta.n_pairs : null)
      if (expected == null && meta2 instanceof Map) expected = meta2.n_pairs
      expected = expected != null ? expected.toString().toInteger() : null

      if (expected == null || expected <= 0) {
          log.warn "[SKIP ERS] ${meta2?.id ?: meta?.id}: missing/invalid n_pairs; will mark incomplete."
          expected = 1
          if (meta2 instanceof Map) meta2.n_pairs = expected
      } else {
          if (meta2 instanceof Map) meta2.n_pairs = expected
      }

      // groupKey: emit as soon as expected BAMS arrive
      def gk = groupKey(meta2.id, expected)

      tuple(gk, meta2, species, ref_name, bam)
  }
  .groupTuple()
  .map { gk, meta2L, speciesL, refL, bams ->        // <-- 解包回单值
      def meta2    = meta2L[0]
      def species  = speciesL[0]
      def ref_name = refL[0]
      tuple(gk, meta2, species, ref_name, bams)
  }

def ch_complete_ers = ch_bam_grouped
  .filter { gk, meta2, species, ref_name, bams ->
      (bams?.size() ?: 0) == (meta2.n_pairs as int)
  }
  .map { gk, meta2, species, ref_name, bams ->
      tuple(meta2, species, ref_name, bams)
  }

def ch_incomplete_ers = ch_bam_grouped
  .filter { gk, meta2, species, ref_name, bams ->
      (bams?.size() ?: 0) < (meta2.n_pairs as int)
  }
  .map { gk, meta2, species, ref_name, bams ->
      def got = bams?.size() ?: 0
      log.warn "[SKIP ERS] ${meta2.id}: only ${got}/${meta2.n_pairs} BAMs succeeded; skip merge+cram+wdl."
      tuple(meta2, species, ref_name, got as int)
  }

// optional: print skipped ERS list
ch_incomplete_ers.view { meta2, species, ref_name, got ->
  "[SKIPPED] ${meta2.id}\t${species}\t${ref_name}\t${got}/${meta2.n_pairs}"
}

// Merge -> CRAM only for complete ERS
def merged_bams    = MERGE_BAMS(ch_complete_ers)
def full_run_crams = BAM_TO_CRAM(merged_bams)


// ============ Path B: Existing CRAMs (skip alignment) ============

def partial_run_crams = branched.partial_run
    .map { meta, species, ref_name ->
        log.info "[PARTIAL SKIP] Found existing CRAM for '${meta.id}'."

        def cram_path = file("${params.outdir}/${meta.id}/alignment/${meta.id}.cram")
        def crai_path = file("${params.outdir}/${meta.id}/alignment/${meta.id}.cram.crai")

        // Guard: only pass through if both files exist
        if (!cram_path.exists() || !crai_path.exists()) {
            log.warn "[PARTIAL SKIP->FALLBACK] Missing published CRAM/CRAI for '${meta.id}', will require full run."
            return null
        }

        tuple(meta, species, ref_name, cram_path, crai_path)
    }
    .filter { it != null }


// ============ Combine both CRAM streams ============

def ch_final_crams = full_run_crams.mix(partial_run_crams)

    // 后续流程
    GENERATE_CRAM_TSV( ch_final_crams )
    GENERATE_WDL_JSON( GENERATE_CRAM_TSV.out.tsv )
    RUN_WDL_VARIANT_CALLING( GENERATE_WDL_JSON.out.json, ch_cromwell_conf )
}

// --- Workflow-level event handlers for reporting ---
workflow.onComplete {
    log.info "Pipeline completed successfully. Output files are in: ${params.outdir}"
}
workflow.onError {
    log.error """
    ----------------------------------------------------
    PIPELINE FAILED
    ----------------------------------------------------
    See '.nextflow.log' for details and the failing task report.
    ----------------------------------------------------
    """
}


// ====================================================================================
//                                    PROCESSES
// ====================================================================================

// Step 1.1: Download SRA and convert to FASTQ
process DOWNLOAD_FASTQ {
    tag { "Download FASTQ for ${meta.id}" }                                  
    label 'down_task'

    errorStrategy 'ignore'
    maxRetries 2

    input:
    tuple val(meta), val(species_name), val(ref_name)

    output:
    tuple val(meta), val(species_name), val(ref_name),path("fastqs"), path("fastqs/.done")

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # ---- tool checks ----
    if ! command -v curl >/dev/null 2>&1; then
      echo "ERROR: curl not found in PATH." >&2
      exit 2
    fi
    if ! command -v md5sum >/dev/null 2>&1; then
      echo "ERROR: md5sum not found in PATH." >&2
      exit 2
    fi

    work_tmp=\$(mktemp -d)
    mkdir -p "\$work_tmp/fastqs" fastqs

    acc="${meta.id}"
    echo "INFO: Resolving runs for accession: ${meta.id}"

    run_list="\$work_tmp/run_ids.txt"
    : > "\$run_list"

    # ---- accession patterns ----
    is_run="^(SRR|ERR|DRR)[0-9]+\$"
    is_experiment="^(SRX|ERX|DRX)[0-9]+\$"
    is_sample="^(SRS|ERS|DRS)[0-9]+\\\$|^(SAMN|SAMEA)[0-9]+\\\$"
    is_study="^(SRP|ERP|DRP)[0-9]+\$"
    is_bioproject="^PRJ(N|E|D)[A-Z][0-9]+\\\$"

    ena_base="https://www.ebi.ac.uk/ena/portal/api"

    # ---- ENA helper with重试 ----
    fetch_ena() {
      local url="\$1"; local out="\$2"
      local tries=3 i=1
      while (( i <= tries )); do
        if curl -fsSL --retry 3 --retry-delay 2 "\$url" -o "\$out"; then
          return 0
        fi
        echo "WARN: ENA request failed (attempt \$i/\$tries): \$url" >&2
        sleep \$((2*i))
        ((i++))
      done
      return 1
    }

    # ---- 根据 accession 类型解析 run 列表 ----
    if [[ "\$acc" =~ \$is_run ]]; then
      printf "%s\\n" "\$acc" > "\$run_list"
    else
      # 先尝试 filereport 直接返回 read_run
      url1="\${ena_base}/filereport?accession=\${acc}&result=read_run&fields=run_accession&format=tsv"
      if fetch_ena "\$url1" "\$work_tmp/filereport.tsv"; then
        tail -n +2 "\$work_tmp/filereport.tsv" | awk -F '\\t' 'NF>0 && \$1!="" {print \$1}' | sort -u >> "\$run_list" || true
      fi

      # 若没拿到 run，再用 search endpoint
      if [[ ! -s "\$run_list" ]]; then
        echo "INFO: filereport returned no runs; trying search endpoint..." >&2
        if   [[ "\$acc" =~ \$is_experiment ]]; then query="experiment_accession%3D%22\${acc}%22"
        elif [[ "\$acc" =~ \$is_sample     ]]; then query="sample_accession%3D%22\${acc}%22"
        elif [[ "\$acc" =~ \$is_study      ]]; then query="study_accession%3D%22\${acc}%22"
        elif [[ "\$acc" =~ \$is_bioproject  ]]; then query="bioproject%3D%22\${acc}%22"
        else
          enc="\${acc}"
          query="(experiment_accession%3D%22\${enc}%22)OR(sample_accession%3D%22\${enc}%22)OR(study_accession%3D%22\${enc}%22)OR(bioproject%3D%22\${enc}%22)"
        fi
        url2="\${ena_base}/search?result=read_run&query=\${query}&fields=run_accession&format=tsv"
        if fetch_ena "\$url2" "\$work_tmp/search.tsv"; then
          tail -n +2 "\$work_tmp/search.tsv" | awk -F '\\t' 'NF>0 && \$1!="" {print \$1}' | sort -u >> "\$run_list" || true
        fi
      fi
    fi

    if [[ ! -s "\$run_list" ]]; then
      echo "ERROR: Could not resolve any SRR/ERR/DRR runs for '${meta.id}' via ENA." >&2
      rm -rf "\$work_tmp"
      exit 1
    fi

    echo -n "INFO: Runs resolved: "
    tr '\\n' ' ' < "\$run_list"; echo

    # ---- 下载 + md5 校验的 helper ----
    download_with_md5() {
      local url="\$1"
      local expected_md5="\$2"
      local out="\$3"

      if [[ -z "\$expected_md5" ]]; then
        echo "ERROR: No MD5 provided for \$url; aborting." >&2
        return 1
      fi

      local tries=3 i=1
      while (( i <= tries )); do
        echo "INFO: [\$i/\$tries] Downloading: \$url"
        # 先下到临时文件，校验通过再 mv
        if ! curl -fSL --retry 3 --retry-delay 2 "\$url" -o "\${out}.tmp"; then
          echo "WARN: Download failed for \$url (attempt \$i/\$tries)" >&2
          ((i++))
          sleep \$((2*i))
          continue
        fi

        local got_md5
        got_md5=\$(md5sum "\${out}.tmp" | awk '{print \$1}')
        if [[ "\$got_md5" == "\$expected_md5" ]]; then
          mv "\${out}.tmp" "\$out"
          echo "INFO: MD5 OK for \$(basename "\$out") (\$got_md5)"
          return 0
        else
          echo "WARN: MD5 mismatch for \$(basename "\$out"): expected=\$expected_md5 got=\$got_md5" >&2
          rm -f "\${out}.tmp"
          ((i++))
          sleep \$((2*i))
        fi
      done

      echo "ERROR: Failed to download with correct MD5 after \$tries attempts: \$url" >&2
      return 1
    }

    # ---- 针对每个 run，通过 filereport 拿 fastq_ftp / fastq_md5，并逐个下载 ----
    while read -r ra; do
      [[ -n "\$ra" ]] || continue
      if [[ ! "\$ra" =~ \$is_run ]]; then
        echo "WARN: Skipping non-run accession in list: \$ra" >&2
        continue
      fi

      echo "INFO: Querying ENA FASTQ URLs for run: \$ra"
      run_report="\$work_tmp/\${ra}_fastq.tsv"
      url_run="\${ena_base}/filereport?accession=\${ra}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5&format=tsv"
      if ! fetch_ena "\$url_run" "\$run_report"; then
        echo "ERROR: Failed to fetch FASTQ report for run \$ra from ENA." >&2
        rm -rf "\$work_tmp"
        exit 2
      fi

      # 解析 fastq_ftp / fastq_md5（分号分隔，多文件对应）
      fq_ftp_col=\$(tail -n +2 "\$run_report" | awk -F '\\t' 'NF>=3 {print \$2; exit}')
      fq_md5_col=\$(tail -n +2 "\$run_report" | awk -F '\\t' 'NF>=3 {print \$3; exit}')

      if [[ -z "\$fq_ftp_col" || -z "\$fq_md5_col" ]]; then
        echo "ERROR: ENA report for run \$ra lacks fastq_ftp or fastq_md5; cannot perform integrity check." >&2
        rm -rf "\$work_tmp"
        exit 3
      fi

      IFS=';' read -r -a fq_urls <<< "\$fq_ftp_col"
      IFS=';' read -r -a fq_md5s <<< "\$fq_md5_col"

      if (( \${#fq_urls[@]} == 0 )); then
        echo "ERROR: No FASTQ URLs found for run \$ra." >&2
        rm -rf "\$work_tmp"
        exit 4
      fi
      if (( \${#fq_urls[@]} != \${#fq_md5s[@]} )); then
        echo "ERROR: Mismatch between number of FASTQ URLs and MD5 entries for run \$ra." >&2
        rm -rf "\$work_tmp"
        exit 5
      fi

      for idx in "\${!fq_urls[@]}"; do
        url="\${fq_urls[\$idx]}"
        md5="\${fq_md5s[\$idx]}"

        # 补全 scheme
        if [[ "\$url" != http*://* ]]; then
          url="https://\${url}"
        fi

        fname=\$(basename "\$url")
        out_path="\$work_tmp/fastqs/\${fname}"

        download_with_md5 "\$url" "\$md5" "\$out_path"
      done

    done < "\$run_list"

    # ---- 将验证过的 fastq.gz 落盘到输出目录，并写一个 MD5SUMS 供后续使用 ----
    echo "INFO: Moving verified FASTQ files to output directory"
    rsync -a "\$work_tmp/fastqs/" fastqs/

    echo "INFO: Writing MD5SUMS file under fastqs/"
    (
      cd fastqs
      # 对输出目录中所有 fastq.gz 再算一遍 md5，方便下游快速复查
      find . -type f -name '*.fastq.gz' -printf '%P\\n' | sort | while read -r f; do
        md5sum "\$f"
      done > MD5SUMS
    )

    touch fastqs/.done
    sync fastqs/.done

    rm -rf "\$work_tmp"
    echo "INFO: DOWNLOAD_FASTQ finished for ${meta.id}"
    """
}

process ALIGN_AND_SORT {
    tag { "Align & Sort ${meta.id} (Pair: ${meta.pair_id})" }                
    label 'alignment_related'

    input:
    tuple val(meta), val(species_name), val(ref_name), path(r1), path(r2)

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}.${meta.pair_id}.sorted.bam"), emit: bam

    script:
    // FIX: 直接在 Bash 内用单引号构造 RG 字符串，避免多重转义
    def ref_file = "${params.global_ref_dir}/${ref_name}.fasta"
    def unsorted_bam = "${meta.id}.${meta.pair_id}.unsorted.bam"            
    def sorted_bam = "${meta.id}.${meta.pair_id}.sorted.bam"
    def sort_mem_per_thread = task.memory ? "${(task.memory.toGiga() * 0.65 / task.cpus).toInteger()}G" : "2G"
    """
    #!/usr/bin/env bash
    set -euo pipefail
    echo "INFO: Aligning reads for ${meta.id}, pair ${meta.pair_id}..."

    # ===== choose fast local tmp =====
    tmp="\$(mktemp -d -p /tmp "\${USER}_nf_\${SLURM_JOB_ID:-nojob}_\${SLURM_ARRAY_TASK_ID:-0}_\${SLURM_PROCID:-0}_XXXXXX")"
    echo "INFO: Using tmp dir: \$tmp"
    df -h "\$tmp" || true

    cleanup() { rm -rf "\$tmp" || true; }
    trap cleanup EXIT

    # ===== alignment + sort =====
    bwa mem -K 100000000 -v 3 -t ${task.cpus} -Y \\
      -R "@RG\\tID:${meta.id}.${meta.pair_id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}" \\
      ${ref_file} ${r1} ${r2} | \\
    samtools sort -@ ${task.cpus} -m ${sort_mem_per_thread} -T "\$tmp/${meta.id}.${meta.pair_id}" -o ${meta.id}.${meta.pair_id}.sorted.bam -

    # 打印一条日志方便在 .command.log 中查看实际分配了多少内存
    echo "INFO: Alignment complete. Sorted with -m ${sort_mem_per_thread} per thread."
    """
}

process MERGE_BAMS {
    tag { "Merge BAMs for ${meta.id}" }
    label 'merge_related'

    input:
    tuple val(meta), val(species_name), val(ref_name), path(bams)                           // FIX: bams 为 list<path>

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}.merged.bam"), path("${meta.id}.merged.bam.bai"), emit: merged_bam

    when:
    bams && bams.size() > 0                                                  // FIX: 无 reads 时直接跳过

    script:
    // FIX: 安全拼接路径列表
    def bam_list = bams.collect{ it.toString() }.join(' ')
    """
    #!/usr/bin/env bash
    set -euo pipefail
    echo "INFO: Merging ${bams.size()} BAM files for ${meta.id}..."
    samtools merge -f -@ ${task.cpus} ${meta.id}.merged.bam ${bam_list}
    echo "INFO: Indexing merged BAM file..."
    samtools index -@ ${task.cpus} ${meta.id}.merged.bam
    """
}

process BAM_TO_CRAM {
    tag { "BAM to CRAM for ${meta.id}" }
    label 'alignment_related'
    publishDir "${params.outdir}/${meta.id}/alignment", mode: 'copy', pattern: "*.{cram,crai}"

    input:
    tuple val(meta), val(species_name), val(ref_name), path(bam), path(bai)

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}.cram"), path("${meta.id}.cram.crai"), emit: cram

    script:
    def ref_file = "${params.global_ref_dir}/${ref_name}.fasta"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    tmp="\$(mktemp -d -p /tmp "\${USER}_nf_\${SLURM_JOB_ID:-nojob}_\${SLURM_ARRAY_TASK_ID:-0}_${meta.id}_cram_XXXXXX")"
    cleanup() { rm -rf "\$tmp" || true; }
    trap cleanup EXIT

    echo "INFO: Converting BAM to CRAM using ${task.cpus} threads; tmp=\$tmp"
    df -h "\$tmp" || true

    samtools view -@ ${task.cpus} -T ${ref_file} -O cram,version=3.0 \\
      -o "\$tmp/${meta.id}.cram" ${bam}

    samtools index -@ ${task.cpus} "\$tmp/${meta.id}.cram"

    mv -f -- "\$tmp/${meta.id}.cram" .
    mv -f -- "\$tmp/${meta.id}.cram.crai" .

    echo "INFO: Deleting intermediate BAM and BAI files..."
    rm -f ${bam} ${bai}
    """
}

process GENERATE_CRAM_TSV {
    tag { "Generate CRAM TSV for ${meta.id}" }
    label 'generation_related'
    publishDir "${params.outdir}/${meta.id}/wdl_inputs", mode: 'copy'

    input:
    tuple val(meta), val(species_name), val(ref_name), path(cram), path(crai)

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}_cram_list.tsv"), emit: tsv

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    #samtools index -@ 2 ${cram}  # 确保有 index
    abspath() {
      if command -v readlink >/dev/null 2>&1; then
        readlink -f "\$1"
      else
        realpath "\$1"
      fi
    }
    cram_path=\$(abspath ${cram})
    crai_path=\$(abspath ${crai})

    echo -e "\${cram_path}\\t\${crai_path}" > "${meta.id}_cram_list.tsv"
    """
}

process GENERATE_WDL_JSON {
    tag { "Generate WDL JSON for ${meta.id}" }
    label 'generation_related'
    publishDir "${params.outdir}/${meta.id}/wdl_inputs", mode: 'copy'

    input:
    tuple val(meta), val(species_name), val(ref_name), path(cram_tsv)

    output:
    tuple val(meta), path("${meta.id}_wdl_inputs.json"), emit: json

    script:
    def mt_fasta = "${params.ref_dir}/${ref_name}.fasta"
    """
    #!/usr/bin/env bash
    set -euo pipefail
    
    ensure_fai() {
        local fa="\$1"; local fai="\${fa}.fai"
        if [[ ! -s "\$fai" ]]; then samtools faidx "\$fa"; fi
        echo "\$fai"
    }
    choose_mt_contig() {
        # Try chrM / MT / or first contig
        local fai="\$1"
        local name=\$(awk '\$1~/^(chr)?M(T)?\$/{print \$1; exit}' "\$fai")
        if [[ -z "\$name" ]]; then
            name=\$(awk 'NR==1{print \$1}' "\$fai")
        fi
        echo "\$name"
    }
    contig_length_from_fai() { awk -v C="\$2" '\$1==C{print \$2}' "\$1"; }

    MT_FAI=\$(ensure_fai "${mt_fasta}")
    CHR_NAME=\$(choose_mt_contig "\$MT_FAI")
    MT_LEN=\$(contig_length_from_fai "\$MT_FAI" "\$CHR_NAME")

    cat > ${meta.id}_wdl_inputs.json <<-EOF
	{
	  "MitochondriaMultiSamplePipeline.picard": "${params.picard_jar}",
	  "MitochondriaMultiSamplePipeline.haplocheckCLI": "${params.haplocheck_jar}",
	  "MitochondriaMultiSamplePipeline.gatk": "${params.gatk_jar}",
	  "MitochondriaMultiSamplePipeline.compress_output_vcf": true,
	  "MitochondriaMultiSamplePipeline.inputSamplesFile": "\$(readlink -f ${cram_tsv})",

	  "MitochondriaMultiSamplePipeline.ref_fasta": "${params.global_ref_dir}/${ref_name}.fasta",
	  "MitochondriaMultiSamplePipeline.ref_fasta_index": "${params.global_ref_dir}/${ref_name}.fasta.fai",
	  "MitochondriaMultiSamplePipeline.ref_dict": "${params.global_ref_dir}/${ref_name}.dict",

	  "MitochondriaMultiSamplePipeline.mt_dict": "${params.ref_dir}/${ref_name}.dict",
	  "MitochondriaMultiSamplePipeline.mt_fasta": "${params.ref_dir}/${ref_name}.fasta",
	  "MitochondriaMultiSamplePipeline.mt_fasta_index": "${params.ref_dir}/${ref_name}.fasta.fai",
	  "MitochondriaMultiSamplePipeline.mt_amb": "${params.ref_dir}/${ref_name}.fasta.amb",
	  "MitochondriaMultiSamplePipeline.mt_ann": "${params.ref_dir}/${ref_name}.fasta.ann",
	  "MitochondriaMultiSamplePipeline.mt_bwt": "${params.ref_dir}/${ref_name}.fasta.bwt",
	  "MitochondriaMultiSamplePipeline.mt_pac": "${params.ref_dir}/${ref_name}.fasta.pac",
	  "MitochondriaMultiSamplePipeline.mt_sa": "${params.ref_dir}/${ref_name}.fasta.sa",

	  "MitochondriaMultiSamplePipeline.mt_shifted_dict": "${params.ref_shift_8000_dir}/${ref_name}.dict",
	  "MitochondriaMultiSamplePipeline.mt_shifted_fasta": "${params.ref_shift_8000_dir}/${ref_name}.fasta",
	  "MitochondriaMultiSamplePipeline.mt_shifted_fasta_index": "${params.ref_shift_8000_dir}/${ref_name}.fasta.fai",
	  "MitochondriaMultiSamplePipeline.mt_shifted_amb": "${params.ref_shift_8000_dir}/${ref_name}.fasta.amb",
	  "MitochondriaMultiSamplePipeline.mt_shifted_ann": "${params.ref_shift_8000_dir}/${ref_name}.fasta.ann",
	  "MitochondriaMultiSamplePipeline.mt_shifted_bwt": "${params.ref_shift_8000_dir}/${ref_name}.fasta.bwt",
	  "MitochondriaMultiSamplePipeline.mt_shifted_pac": "${params.ref_shift_8000_dir}/${ref_name}.fasta.pac",
	  "MitochondriaMultiSamplePipeline.mt_shifted_sa": "${params.ref_shift_8000_dir}/${ref_name}.fasta.sa",

	  "MitochondriaMultiSamplePipeline.shift_back_chain": "${params.shift_back_chain_dir}/${ref_name}_ShiftBack.chain",
	  "MitochondriaMultiSamplePipeline.non_control_region_interval_list": "${params.ref_interval_dir}/${ref_name}_non_control_region.interval_list",
	  "MitochondriaMultiSamplePipeline.control_region_shifted_reference_interval_list": "${params.ref_interval_dir}/${ref_name}_control_region_shifted.interval_list",

	  "MitochondriaMultiSamplePipeline.mt_chr_name": "\$CHR_NAME",
	  "MitochondriaMultiSamplePipeline.mt_length": \$MT_LEN,
	  "MitochondriaMultiSamplePipeline.mt_nc_start": ${params.mt_nc_start},
	  "MitochondriaMultiSamplePipeline.mt_right_pad": ${params.mt_right_pad},
	  "MitochondriaMultiSamplePipeline.mt_shift": ${params.mt_shift}
	}
	EOF
    """
}

process RUN_WDL_VARIANT_CALLING {
    tag "Variant Calling on ${meta.id}"
    label 'wdl_related'
    // Also update the publishDir pattern to copy all VCF-related files
    publishDir "${params.outdir}/${meta.id}/variant_calling/", mode: 'copy'

    input:
    tuple val(meta), val(wdl_inputs_json)
    path cromwell_config

    output:
    path "cromwell-executions/**"

    script:
    """
    #!/bin/bash
    set -e

    # ... (The cat commands for cromwell_options.json and sbatch_throttle.sh are correct and remain unchanged) ...
    cat > cromwell_options.json <<EOF
    {
      "final_workflow_outputs_dir": ".",
      "final_workflow_log_dir": ".",
      "default_runtime_attributes": {
        "queue": "${params.cromwell_options.queue ?: ''}",
        "cpus": ${params.cromwell_options.cpus},
        "memory": ${params.cromwell_options.memory},
        "runtime_minutes": ${params.cromwell_options.runtime_minutes}
      }
    }
    EOF
    cat > sbatch_throttle.sh <<'EOS'
    #!/usr/bin/env bash
    set -euo pipefail
    PER_HOUR="\${SUBMITS_PER_HOUR:-180}"
    (( PER_HOUR > 0 )) || PER_HOUR=180
    MIN_GAP=\$(( 3600 / PER_HOUR ))
    STATE_DIR="\${HOME}/.sbatch_rate"
    LOCK_FILE="\${STATE_DIR}/lock"
    TS_FILE="\${STATE_DIR}/last_submit.ts"
    mkdir -p "\${STATE_DIR}"
    exec 200>"\${LOCK_FILE}"
    flock 200
    now=\$(date +%s); last=0
    [[ -f "\${TS_FILE}" ]] && read -r last < "\${TS_FILE}" || true
    delta=\$(( now - last ))
    if (( delta < MIN_GAP )); then
      sleep \$(( MIN_GAP - delta ))
    fi
    date +%s > "\${TS_FILE}"
    unset SLURM_CONF || true
    exec sbatch "\$@"
    EOS
    chmod +x sbatch_throttle.sh

    export SUBMITS_PER_HOUR="${params.cromwell_submit_rate_limit ?: '180'}"

    echo "INFO: Forcing execution with Java from JAVA_HOME: \$JAVA_HOME"
    "\${JAVA_HOME}/bin/java" -Xmx4G -Dconfig.file=${cromwell_config} \\
         -jar ${params.cromwell_jar} run \\
         ${params.wdl_script} \\
         --inputs ${wdl_inputs_json} \\
         --options cromwell_options.json
    
    """
}
