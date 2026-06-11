#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * PRIMATE mtDNA PIPELINE - NUMT DECOY ROUND1
 * ------------------------------------------------------------
 * Input sample TSV columns:
 *   col1 = sample_id
 *   col2 = species_name
 *   col3 = ref_name
 *
 * This pipeline starts from existing WGS CRAMs and performs:
 *   1) locate CRAM/CRAI
 *   2) build chrM + NUMT decoy reference
 *   3) MitoHPC-like extraction of primary reads overlapping chrM + NUMT intervals
 *      and retain only templates with both mates in the target subset
 *   4) FASTQ conversion + realignment to chrM+NUMT decoy
 *   6) keep only final chrM alignments
 *   7) convert cleaned chrM BAM -> CRAM
 *   8) feed into existing RUN_WDL_VARIANT_CALLING module
 *
 * Notes:
 *   - NUMT BED is treated as decoy only.
 *   - No NUMT variant calling is done here.
 *   - Output CRAM contains chrM-only cleaned reads aligned to chrM.
 */

log.info """
PRIMATE NUMT-DECOY ROUND1 PIPELINE START
========================================
Sample TSV:              ${params.sample_tsv}
Output Directory:        ${params.outdir}
CRAM search dirs:        ${params.cram_dirs}
Whole-genome ref dir:    ${params.global_ref_dir}
mt ref dir:              ${params.ref_dir}
NUMT BED dir:            ${params.numt_bed_dir}
WDL Script:              ${params.wdl_script}
========================================
"""

ch_samples = Channel.fromPath(params.sample_tsv)
    .splitCsv(header: false, sep: '\t')
    .filter { row -> row.size() >= 3 && row[0]?.trim() && row[1]?.trim() && row[2]?.trim() }
    .map { row ->
        def sample_id = row[0].trim()
        def species = row[1].trim()
        def ref_name = row[2].trim()

        def round1_dir = file("${params.outdir}/${sample_id}/round_1")
        def vc_dir     = file("${params.outdir}/${sample_id}/round_1_variant_calling_decoy")

        if (round1_dir.exists() && vc_dir.exists()) {
            log.info "SKIP completed sample ${sample_id}: ${round1_dir} and ${vc_dir} already exist"
            return null
        }

        def meta = [id: sample_id]
        tuple(meta, species, ref_name)
    }
    .filter { it != null }

ch_cromwell_conf = file("${baseDir}/cromwell.conf")

workflow {

    LOCATE_CRAM(ch_samples)
    PREPARE_DECOY_REFERENCE(LOCATE_CRAM.out.cram_info)

    // MitoHPC-like candidate read selection:
    //   1) subset primary alignments overlapping chrM + NUMT intervals
    //   2) name-sort
    //   3) keep only templates with both mates present in the target subset
    // This avoids rescuing arbitrary NUMT/chrM mates from elsewhere in the genome.
    EXTRACT_CANDIDATE_READS(LOCATE_CRAM.out.cram_info)

    ch_bam_plus_ref = EXTRACT_CANDIDATE_READS.out.bam_with_mates.join(PREPARE_DECOY_REFERENCE.out.decoy_ref, by: [0,1,2])
    REALIGN_TO_DECOY(ch_bam_plus_ref)

    EXTRACT_FINAL_CHRM_BAM(REALIGN_TO_DECOY.out.realign_bam)

    ch_final_bam_plus_cram = EXTRACT_FINAL_CHRM_BAM.out.final_bam.join(LOCATE_CRAM.out.cram_info, by: [0,1,2])
    BAM_TO_CRAM_CLEAN(ch_final_bam_plus_cram)

    GENERATE_CRAM_TSV(BAM_TO_CRAM_CLEAN.out.cram)
    GENERATE_WDL_JSON(GENERATE_CRAM_TSV.out.tsv)
    RUN_WDL_VARIANT_CALLING(GENERATE_WDL_JSON.out.json, ch_cromwell_conf)
}

workflow.onComplete {
    if( workflow.success ) {
        log.info "Pipeline completed successfully. Output files are in: ${params.outdir}"
    } else {
        log.error "Pipeline finished with errors. Check .nextflow.log and failed task work dirs."
    }
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

process LOCATE_CRAM {
    tag { "Locate CRAM for ${meta.id}" }
    label 'generation_related'

    input:
    tuple val(meta), val(species_name), val(ref_name)

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}.input.cram"), path("${meta.id}.input.cram.crai"), emit: cram_info

    script:
    def cramDirs = (params.cram_dirs instanceof List ? params.cram_dirs : params.cram_dirs.toString().split(',')*.trim()).findAll{ it }
    def cramDirsBash = cramDirs.collect { '"' + it + '"' }.join(' ')
    """
    #!/usr/bin/env bash
    set -euo pipefail

    sample_id="${meta.id}"
    found_cram=""
    found_crai=""

    for d in ${cramDirsBash}; do
      cand_cram="\${d}/\${sample_id}/alignment/\${sample_id}.cram"
      cand_crai="\${d}/\${sample_id}/alignment/\${sample_id}.cram.crai"
      if [[ -s "\${cand_cram}" && -s "\${cand_crai}" ]]; then
        found_cram="\${cand_cram}"
        found_crai="\${cand_crai}"
        break
      fi
    done

    if [[ -z "\${found_cram}" ]]; then
      echo "ERROR: CRAM/CRAI not found for \${sample_id} in any search dir." >&2
      exit 1
    fi

    ln -sf "\${found_cram}" ${meta.id}.input.cram
    ln -sf "\${found_crai}" ${meta.id}.input.cram.crai
    """
}

process PREPARE_DECOY_REFERENCE {
    tag { "Prepare decoy ref for ${meta.id}" }
    label 'alignment_related'
    publishDir "${params.outdir}/${meta.id}/round_1/numt_decoy_ref", mode: 'copy'

    input:
    tuple val(meta), val(species_name), val(ref_name), path(cram), path(crai)

    output:
    // Keep the original decoy_ref output structure unchanged,
    // so downstream processes using decoy_ref will not break.
    tuple val(meta), val(species_name), val(ref_name),
          path("${meta.id}.chrM_plus_numt.fa"),
          path("${meta.id}.chrM_plus_numt.fa.fai"),
          path("${meta.id}.chrM_plus_numt.fa.amb"),
          path("${meta.id}.chrM_plus_numt.fa.ann"),
          path("${meta.id}.chrM_plus_numt.fa.bwt"),
          path("${meta.id}.chrM_plus_numt.fa.pac"),
          path("${meta.id}.chrM_plus_numt.fa.sa"),
          emit: decoy_ref

    // New output for Round 2:
    // original NUMT sequences extracted from the original whole-genome reference.
    tuple val(meta), val(species_name), val(ref_name),
          path("${meta.id}.original_numt.fa"),
          path("${meta.id}.original_numt.fa.fai"),
          emit: original_numt_fa

    script:
    def whole_ref = "${params.global_ref_dir}/${ref_name}.fasta"
    def numt_bed  = "${params.numt_bed_dir}/${meta.id}${params.numt_bed_suffix}"
    def mt_contig = params.mt_contig ?: "chrM"

    """
    #!/usr/bin/env bash
    set -euo pipefail

    REF="${whole_ref}"
    BED="${numt_bed}"
    MT_CONTIG="${mt_contig}"

    echo "[INFO] Sample: ${meta.id}"
    echo "[INFO] Species: ${species_name}"
    echo "[INFO] Reference name: ${ref_name}"
    echo "[INFO] Whole-genome reference: \${REF}"
    echo "[INFO] NUMT BED: \${BED}"
    echo "[INFO] mtDNA contig: \${MT_CONTIG}"

    [[ -s "\${REF}" ]] || { echo "ERROR: Missing reference fasta: \${REF}" >&2; exit 1; }
    [[ -s "\${REF}.fai" ]] || { echo "ERROR: Missing reference fasta index: \${REF}.fai" >&2; exit 1; }
    [[ -e "\${BED}" ]] || { echo "ERROR: Missing NUMT BED file: \${BED}" >&2; exit 1; }
    if [[ ! -s "\${BED}" ]]; then
        echo "[WARN] NUMT BED exists but is empty: \${BED}"
        echo "[WARN] Proceeding with chrM-only decoy reference for this sample."
    fi

    # 1. Extract original chrM sequence from whole-genome reference.
    samtools faidx "\${REF}" "\${MT_CONTIG}" > chrM.fa

    [[ -s chrM.fa ]] || {
        echo "ERROR: Failed to extract \${MT_CONTIG} from \${REF}" >&2
        exit 1
    }

    # 2. Prepare named NUMT BED.
    # The name will become the FASTA header when using bedtools getfasta -nameOnly.
    awk 'BEGIN{OFS="\\t"}
         \$0 !~ /^#/ && NF >= 3 && \$3 > \$2 {
           name=sprintf("NUMT_%s_%s_%s_%s", NR, \$1, \$2+1, \$3)
           print \$1,\$2,\$3,name
         }' "\${BED}" > numt.named.bed

    if [[ -s numt.named.bed ]]; then
        # 3. Extract original NUMT sequences from the whole-genome reference.
        bedtools getfasta \
            -fi "\${REF}" \
            -bed numt.named.bed \
            -nameOnly \
            -fo numt_decoy.fa

        [[ -s numt_decoy.fa ]] || {
            echo "ERROR: bedtools getfasta produced empty NUMT fasta despite non-empty parsed BED" >&2
            exit 1
        }
    else
        echo "[WARN] No valid NUMT intervals after parsing BED: \${BED}"
        echo "[WARN] Creating an empty original NUMT FASTA and chrM-only decoy reference."
        : > numt_decoy.fa
    fi

    # 4. Publish original NUMT FASTA for Round 2 mtSwirl-like reference:
    #    consensus chrM + original NUMT. If BED is empty, this FASTA is intentionally empty.
    cp numt_decoy.fa ${meta.id}.original_numt.fa
    if [[ -s ${meta.id}.original_numt.fa ]]; then
        samtools faidx ${meta.id}.original_numt.fa
    else
        : > ${meta.id}.original_numt.fa.fai
    fi

    echo "[INFO] Number of original NUMT contigs:"
    grep -c '^>' ${meta.id}.original_numt.fa || true

    # 5. Build Round 1 decoy reference: original chrM + original NUMT.
    #    If no valid NUMT intervals exist, numt_decoy.fa is empty and this becomes chrM-only.
    cat chrM.fa numt_decoy.fa > ${meta.id}.chrM_plus_numt.fa

    [[ -s ${meta.id}.chrM_plus_numt.fa ]] || {
        echo "ERROR: Failed to create ${meta.id}.chrM_plus_numt.fa" >&2
        exit 1
    }

    # 6. Index decoy reference.
    samtools faidx ${meta.id}.chrM_plus_numt.fa
    bwa index ${meta.id}.chrM_plus_numt.fa

    echo "[INFO] Decoy reference created successfully:"
    ls -lh ${meta.id}.chrM_plus_numt.fa*
    echo "[INFO] Original NUMT FASTA created successfully:"
    ls -lh ${meta.id}.original_numt.fa*
    """
}

process EXTRACT_CANDIDATE_READS {
    tag { "MitoHPC-like target-paired reads for ${meta.id}" }
    label 'alignment_related'
    publishDir "${params.outdir}/${meta.id}/round_1/candidate_reads", mode: 'copy', pattern: "*.{bam,bai,tsv,bed}"

    input:
    tuple val(meta), val(species_name), val(ref_name), path(cram), path(crai)

    output:
    // This BAM replaces the previous chrM/NUMT seed + arbitrary mate-rescue BAM.
    // It contains only primary records where both mates were present in the
    // chrM + NUMT target subset, following the MitoHPC-style restriction.
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}.with_mates.bam"), emit: bam_with_mates

    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}.with_mates.bam.bai"), emit: bam_with_mates_index

    tuple val(meta), val(species_name), val(ref_name),
          path("${meta.id}.target_chrM_plus_NUMT.bed"),
          path("${meta.id}.candidate_read_selection.summary.tsv"),
          emit: candidate_summary

    script:
    def mt_contig = params.mt_contig ?: "chrM"
    def numt_bed  = "${params.numt_bed_dir}/${meta.id}${params.numt_bed_suffix}"

    """
    #!/usr/bin/env bash
    set -euo pipefail

    SAMPLE_ID="${meta.id}"
    REF="${params.global_ref_dir}/${ref_name}.fasta"
    BED="${numt_bed}"
    MT_CONTIG="${mt_contig}"
    THREADS=${task.cpus ?: 4}

    echo "[INFO] Sample: \${SAMPLE_ID}"
    echo "[INFO] MitoHPC-like candidate read extraction"
    echo "[INFO] Strategy: primary reads overlapping chrM + NUMT intervals, then keep only read-pairs with both mates in target subset"
    echo "[INFO] Reference: \${REF}"
    echo "[INFO] NUMT BED: \${BED}"
    echo "[INFO] mtDNA contig: \${MT_CONTIG}"

    [[ -s "\${REF}" ]] || { echo "ERROR: Missing whole-genome reference: \${REF}" >&2; exit 1; }
    [[ -s "\${REF}.fai" ]] || { echo "ERROR: Missing whole-genome reference index: \${REF}.fai" >&2; exit 1; }
    [[ -s "${cram}" ]] || { echo "ERROR: Missing input CRAM: ${cram}" >&2; exit 1; }
    [[ -e "\${BED}" ]] || { echo "ERROR: Missing NUMT BED file: \${BED}" >&2; exit 1; }
    if [[ ! -s "\${BED}" ]]; then
        echo "[WARN] NUMT BED exists but is empty: \${BED}"
        echo "[WARN] Candidate target BED will contain chrM only."
    fi

    MT_LEN=\$(awk -v C="\${MT_CONTIG}" '\$1==C{print \$2; exit}' "\${REF}.fai")
    [[ -n "\${MT_LEN}" ]] || { echo "ERROR: Cannot find \${MT_CONTIG} in \${REF}.fai" >&2; exit 1; }

    # Build target BED: full chrM plus sample-specific NUMT intervals.
    # BED is 0-based half-open. chrM is added as 0..MT_LEN.
    awk -v mt="\${MT_CONTIG}" -v len="\${MT_LEN}" 'BEGIN{OFS="\\t"; print mt,0,len}' > \${SAMPLE_ID}.target_chrM_plus_NUMT.bed
    awk 'BEGIN{OFS="\\t"} \$0 !~ /^#/ && NF>=3 && \$3>\$2 {print \$1,\$2,\$3}' "\${BED}" >> \${SAMPLE_ID}.target_chrM_plus_NUMT.bed

    [[ -s \${SAMPLE_ID}.target_chrM_plus_NUMT.bed ]] || { echo "ERROR: target BED is empty" >&2; exit 1; }

    echo "[INFO] Target intervals:"
    wc -l \${SAMPLE_ID}.target_chrM_plus_NUMT.bed

    # Step 1: subset primary alignments overlapping chrM + NUMT intervals.
    # -F 0x900 matches the MitoHPC-style exclusion of secondary and supplementary records.
    samtools view \
        -@ "\${THREADS}" \
        -T "\${REF}" \
        -b \
        -F 0x900 \
        -L \${SAMPLE_ID}.target_chrM_plus_NUMT.bed \
        "${cram}" \
        > \${SAMPLE_ID}.target.primary.unsorted.bam

    target_records=\$(samtools view -c \${SAMPLE_ID}.target.primary.unsorted.bam || echo 0)
    echo "[INFO] Primary records overlapping chrM+NUMT target intervals: \${target_records}"
    [[ "\${target_records}" -gt 0 ]] || { echo "ERROR: no primary target records extracted" >&2; exit 1; }

    # Step 2: name-sort, then keep only read names with both mates present in the target subset.
    # This is the key difference from the previous arbitrary mate-rescue strategy.
    samtools sort -n -@ "\${THREADS}" \
        -o \${SAMPLE_ID}.target.primary.name.bam \
        \${SAMPLE_ID}.target.primary.unsorted.bam

    samtools view -h \${SAMPLE_ID}.target.primary.name.bam | \
      awk 'BEGIN{FS=OFS="\\t"}
           /^@/ {print; next}
           {
             q=\$1
             if (curr=="") {curr=q; n=0}
             if (q!=curr) {
               if (n>=2) {for(i=1;i<=n;i++) print rec[i]}
               delete rec; n=0; curr=q
             }
             rec[++n]=\$0
           }
           END{
             if (n>=2) {for(i=1;i<=n;i++) print rec[i]}
           }' | \
      samtools view -@ "\${THREADS}" -b -o \${SAMPLE_ID}.target_pairs.name.bam -

    pair_records=\$(samtools view -c \${SAMPLE_ID}.target_pairs.name.bam || echo 0)
    echo "[INFO] Records retained after requiring both mates in chrM+NUMT target subset: \${pair_records}"
    [[ "\${pair_records}" -gt 0 ]] || { echo "ERROR: no complete target-region read pairs retained" >&2; exit 1; }

    # Coordinate-sort for downstream handling. REALIGN_TO_DECOY will name-sort again before FASTQ conversion.
    samtools sort -@ "\${THREADS}" \
        -o \${SAMPLE_ID}.with_mates.bam \
        \${SAMPLE_ID}.target_pairs.name.bam
    samtools index -@ "\${THREADS}" \${SAMPLE_ID}.with_mates.bam

    # Diagnostic summary.
    target_qnames=\$(samtools view \${SAMPLE_ID}.target.primary.unsorted.bam | cut -f1 | sort -u | wc -l)
    pair_qnames=\$(samtools view \${SAMPLE_ID}.with_mates.bam | cut -f1 | sort -u | wc -l)
    chrm_records=\$(samtools view \${SAMPLE_ID}.with_mates.bam | awk -v mt="\${MT_CONTIG}" '\$3==mt{n++} END{print n+0}')
    numt_records=\$(samtools view \${SAMPLE_ID}.with_mates.bam | awk -v mt="\${MT_CONTIG}" '\$3!=mt{n++} END{print n+0}')

    cat > \${SAMPLE_ID}.candidate_read_selection.summary.tsv <<EOF_SUMMARY
metric\tvalue
target_intervals_n\t\$(wc -l < \${SAMPLE_ID}.target_chrM_plus_NUMT.bed)
primary_target_records_before_pair_filter\t\${target_records}
primary_target_readnames_before_pair_filter\t\${target_qnames}
records_after_requiring_both_mates_in_target_subset\t\${pair_records}
readnames_after_requiring_both_mates_in_target_subset\t\${pair_qnames}
retained_records_mapped_to_chrM_in_original_cram\t\${chrm_records}
retained_records_mapped_to_NUMT_intervals_or_other_target_contigs_in_original_cram\t\${numt_records}
selection_mode\tmitoHPC_like_chrM_plus_NUMT_target_pairs_no_arbitrary_mate_rescue
excluded_flags_initial_subset\t0x900_secondary_and_supplementary
EOF_SUMMARY

    echo "[INFO] Candidate read selection summary:"
    cat \${SAMPLE_ID}.candidate_read_selection.summary.tsv

    rm -f \${SAMPLE_ID}.target.primary.unsorted.bam \${SAMPLE_ID}.target.primary.name.bam \${SAMPLE_ID}.target_pairs.name.bam
    """
}


process REALIGN_TO_DECOY {
    tag { "Realign to decoy for ${meta.id}" }
    label 'alignment_related'
    publishDir "${params.outdir}/${meta.id}/round_1/decoy_realign", mode: 'copy', pattern: "*.decoy_realign.bam*"

    input:
    tuple val(meta), val(species_name), val(ref_name), path(bam_with_mates), path(decoy_fa), path(decoy_fai), path(decoy_amb), path(decoy_ann), path(decoy_bwt), path(decoy_pac), path(decoy_sa)

    output:
    tuple val(meta), val(species_name), val(ref_name),
          path("${meta.id}.decoy_realign.bam"),
          path("${meta.id}.decoy_realign.bam.bai"),
          emit: realign_bam

    script:
    def sort_mem_per_thread = task.memory ? "${(task.memory.toGiga() * 0.65 / Math.max(task.cpus ?: 1,1)).toInteger()}G" : "2G"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    tmp="\$(mktemp -d -p /tmp "\${USER}_nf_${meta.id}_decoy_XXXXXX")"
    cleanup() { rm -rf "\$tmp" || true; }
    trap cleanup EXIT

    samtools sort -n -@ ${task.cpus ?: 4} -o ${meta.id}.with_mates.name.bam ${bam_with_mates}

    samtools fastq \
      -1 ${meta.id}.R1.fastq.gz \
      -2 ${meta.id}.R2.fastq.gz \
      -0 /dev/null -s /dev/null -n \
      ${meta.id}.with_mates.name.bam

    bwa mem -K 100000000 -v 3 -t ${task.cpus ?: 4} \
      -R '@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:${meta.id}\\tPL:ILLUMINA\\tPU:${meta.id}' \
      ${decoy_fa} \
      ${meta.id}.R1.fastq.gz ${meta.id}.R2.fastq.gz | \
      samtools sort -@ ${task.cpus ?: 4} -m ${sort_mem_per_thread} -T "\$tmp/${meta.id}" -o ${meta.id}.decoy_realign.bam -

    samtools index -@ ${task.cpus ?: 4} ${meta.id}.decoy_realign.bam

    rm -f ${meta.id}.with_mates.name.bam ${meta.id}.R1.fastq.gz ${meta.id}.R2.fastq.gz
    """
}

process EXTRACT_FINAL_CHRM_BAM {
    tag { "Extract final primary chrM BAM for ${meta.id}" }
    label 'alignment_related'
    publishDir "${params.outdir}/${meta.id}/round_1/chrM_clean", mode: 'copy', pattern: "*.final_chrM.sorted.bam*"

    input:
    tuple val(meta), val(species_name), val(ref_name),
          path(realign_bam), path(realign_bai)

    output:
    tuple val(meta), val(species_name), val(ref_name),
          path("${meta.id}.final_chrM.sorted.bam"),
          path("${meta.id}.final_chrM.sorted.bam.bai"),
          emit: final_bam

    script:
    def sort_mem_per_thread = task.memory ? "${(task.memory.toGiga() * 0.65 / Math.max(task.cpus ?: 1,1)).toInteger()}G" : "2G"
    def mt_contig = params.mt_contig ?: "chrM"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    MT_CONTIG="${mt_contig}"

    tmp="\$(mktemp -d -p /tmp "\${USER}_nf_${meta.id}_finalchrm_XXXXXX")"
    cleanup() { rm -rf "\$tmp" || true; }
    trap cleanup EXIT

    # Diagnostics before final extraction.
    all_chrM_records=\$(samtools view -c ${realign_bam} "\${MT_CONTIG}" || echo 0)
    primary_chrM_records=\$(samtools view -c -F 0x900 ${realign_bam} "\${MT_CONTIG}" || echo 0)
    secondary_chrM_records=\$(samtools view -c -f 0x100 ${realign_bam} "\${MT_CONTIG}" || echo 0)
    supplementary_chrM_records=\$(samtools view -c -f 0x800 ${realign_bam} "\${MT_CONTIG}" || echo 0)

    echo "[INFO] chrM records in decoy-realigned BAM:"
    echo "[INFO]   all chrM records:           \${all_chrM_records}"
    echo "[INFO]   primary chrM records:       \${primary_chrM_records}"
    echo "[INFO]   secondary chrM records:     \${secondary_chrM_records}"
    echo "[INFO]   supplementary chrM records: \${supplementary_chrM_records}"

    # 1) chrM-only header: keep chrM @SQ, drop NUMT_* @SQ, keep all non-@SQ header lines.
    samtools view -H ${realign_bam} | \
      awk -v mt="\${MT_CONTIG}" 'BEGIN{OFS="\\t"}
           /^@SQ/ {
             if (\$2=="SN:" mt) print;
             next
           }
           { print }' > "\$tmp/chrM_only.header.sam"

    # 2) chrM-only body, primary alignments only.
    # -F 0x900 removes secondary and supplementary records. This prevents reads whose
    # primary alignment is NUMT but secondary/supplementary alignment is chrM from entering calling.
    samtools view -@ ${task.cpus ?: 4} -F 0x900 ${realign_bam} "\${MT_CONTIG}" > "\$tmp/chrM_only.body.sam"

    # 3) rebuild BAM with cleaned chrM-only header.
    cat "\$tmp/chrM_only.header.sam" "\$tmp/chrM_only.body.sam" | \
      samtools view -@ ${task.cpus ?: 4} -b -o "\$tmp/${meta.id}.final_chrM.bam" -

    # 4) sort + index.
    samtools sort -@ ${task.cpus ?: 4} -m ${sort_mem_per_thread} \
      -T "\$tmp/${meta.id}" \
      -o ${meta.id}.final_chrM.sorted.bam \
      "\$tmp/${meta.id}.final_chrM.bam"

    samtools index -@ ${task.cpus ?: 4} ${meta.id}.final_chrM.sorted.bam

    # 5) sanity check.
    samtools idxstats ${meta.id}.final_chrM.sorted.bam > ${meta.id}.final_chrM.sorted.idxstats.txt
    samtools quickcheck -v ${meta.id}.final_chrM.sorted.bam
    """
}

process BAM_TO_CRAM_CLEAN {
    tag { "BAM to cleaned CRAM for ${meta.id}" }
    label 'alignment_related'
    publishDir "${params.outdir}/${meta.id}/round_1/alignment_numt_decoy", mode: 'copy', pattern: "*.{cram,crai}"

    input:
    tuple val(meta), val(species_name), val(ref_name), path(final_bam), path(final_bai), path(orig_cram), path(orig_crai)

    output:
    tuple val(meta), val(species_name), val(ref_name),
          path("${meta.id}.numt_decoy.clean.cram"),
          path("${meta.id}.numt_decoy.clean.cram.crai"),
          emit: cram

    script:
    def mt_ref = "${params.ref_dir}/${ref_name}.fasta"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    tmp="\$(mktemp -d -p /tmp "\${USER}_nf_${meta.id}_cleancram_XXXXXX")"
    cleanup() { rm -rf "\$tmp" || true; }
    trap cleanup EXIT

    samtools view -@ ${task.cpus ?: 4} \
      -T ${mt_ref} \
      -O cram,version=3.0 \
      -o "\$tmp/${meta.id}.numt_decoy.clean.cram" \
      ${final_bam}

    samtools index -@ ${task.cpus ?: 4} "\$tmp/${meta.id}.numt_decoy.clean.cram"

    samtools quickcheck -v "\$tmp/${meta.id}.numt_decoy.clean.cram"
    samtools idxstats "\$tmp/${meta.id}.numt_decoy.clean.cram" > "\$tmp/${meta.id}.numt_decoy.clean.idxstats.txt"

    mv -f "\$tmp/${meta.id}.numt_decoy.clean.cram" .
    mv -f "\$tmp/${meta.id}.numt_decoy.clean.cram.crai" .
    """
}

process GENERATE_CRAM_TSV {
    tag { "Generate CRAM TSV for ${meta.id}" }
    label 'generation_related'
    publishDir "${params.outdir}/${meta.id}/round_1/wdl_inputs", mode: 'copy'

    input:
    tuple val(meta), val(species_name), val(ref_name), path(cram), path(crai)

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}_cram_list.tsv"), emit: tsv

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    abspath() {
      if command -v readlink >/dev/null 2>&1; then
        readlink -f "\$1"
      else
        realpath "\$1"
      fi
    }
    cram_path=\$(abspath ${cram})
    crai_path=\$(abspath ${crai})

    echo -e "\${cram_path}\t\${crai_path}" > "${meta.id}_cram_list.tsv"
    """
}

process GENERATE_WDL_JSON {
    tag { "Generate WDL JSON for ${meta.id}" }
    label 'generation_related'
    publishDir "${params.outdir}/${meta.id}/round_1/wdl_inputs", mode: 'copy'

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

    cat > ${meta.id}_wdl_inputs.json <<-EOFJSON
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
  "MitochondriaMultiSamplePipeline.nuc_interval_list": "${params.numt_bed_dir}/${meta.id}${params.numt_bed_suffix}",
  "MitochondriaMultiSamplePipeline.use_haplotype_caller_nucdna": true,
  "MitochondriaMultiSamplePipeline.haplotype_caller_nucdna_dp_lower_bound": 10,

  "MitochondriaMultiSamplePipeline.mt_chr_name": "\$CHR_NAME",
  "MitochondriaMultiSamplePipeline.mt_length": \$MT_LEN,
  "MitochondriaMultiSamplePipeline.mt_nc_start": ${params.mt_nc_start},
  "MitochondriaMultiSamplePipeline.mt_right_pad": ${params.mt_right_pad},
  "MitochondriaMultiSamplePipeline.mt_shift": ${params.mt_shift}
}
EOFJSON
    """
}

process RUN_WDL_VARIANT_CALLING {
    tag "Variant Calling on ${meta.id}"
    label 'wdl_related'
    publishDir "${params.outdir}/${meta.id}/round_1_variant_calling_decoy/", mode: 'copy'

    input:
    tuple val(meta), path(wdl_inputs_json)
    path cromwell_config

    output:
    path "cromwell-executions/**"

    script:
    """
    #!/bin/bash
    set -e

    cat > cromwell_options.json <<EOFOPT
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
EOFOPT

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
if (( now - last < MIN_GAP )); then
  sleep \$(( MIN_GAP - (now - last) ))
fi
date +%s > "\${TS_FILE}"
unset SLURM_CONF || true
exec sbatch "\$@"
EOS
    chmod +x sbatch_throttle.sh

    export SUBMITS_PER_HOUR="${params.cromwell_submit_rate_limit ?: '180'}"

    "\${JAVA_HOME}/bin/java" -Dconfig.file=${cromwell_config} \
         -jar ${params.cromwell_jar} run \
         ${params.wdl_script} \
         --inputs ${wdl_inputs_json} \
         --options cromwell_options.json
    """
}