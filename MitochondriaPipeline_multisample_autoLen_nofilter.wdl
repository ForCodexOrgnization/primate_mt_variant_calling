version 1.0

import "AlignAndCall.wdl" as AlignAndCall

workflow MitochondriaMultiSamplePipeline {

  meta {
    description: "Takes in primate BAM/CRAM and outputs mitochondrial SNP/INDEL VCFs."
    allowNestedInputs: true
  }

  input {
    # TSV with columns: input_bam  input_bai
    File inputSamplesFile

    # Species-specific / reference-specific parameters (now supplied via JSON)
    String mt_chr_name
    Int    mt_length
    Int    mt_nc_start = 576
    Int    mt_right_pad = 545
    Int    mt_shift = 8000

    # Coverage used by optional NuMT filter inside AlignAndCall
    Float autosomal_coverage = 0

    # Tools / behavior
    String gatk
    String picard
    String haplocheckCLI
    Boolean compress_output_vcf

    # Full reference only required if starting with CRAM
    File? ref_fasta
    File? ref_fasta_index
    File? ref_dict

    # Unshifted mt reference (and indexes)
    File mt_dict
    File mt_fasta
    File mt_fasta_index
    File mt_amb
    File mt_ann
    File mt_bwt
    File mt_pac
    File mt_sa

    # Shifted mt reference (and indexes) for control region
    File mt_shifted_dict
    File mt_shifted_fasta
    File mt_shifted_fasta_index
    File mt_shifted_amb
    File mt_shifted_ann
    File mt_shifted_bwt
    File mt_shifted_pac
    File mt_shifted_sa

    # Liftover chain (shifted -> original)
    File shift_back_chain

    # Picard interval_lists (non-control on unshifted, control on shifted)
    File non_control_region_interval_list
    File control_region_shifted_reference_interval_list
  }

  Array[Array[String]] inputSamples = read_tsv(inputSamplesFile)

  scatter (sample in inputSamples) {

    call SubsetBamToChrM {
      input:
        gatk            = gatk,
        input_bam       = sample[0],
        input_bai       = sample[1],
        contig_name     = mt_chr_name,
        ref_fasta       = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict        = ref_dict
    }

    call RevertSam {
      input:
        picard      = picard,
        input_bam   = SubsetBamToChrM.output_bam
    }

    String base_name = basename(SubsetBamToChrM.output_bam, ".bam")

    call AlignAndCall.AlignAndCall as AlignAndCall {
      input:
        unmapped_bam            = RevertSam.unmapped_bam,
        autosomal_coverage      = autosomal_coverage,
        base_name               = base_name,
        picard                  = picard,
        gatk                    = gatk,
        haplocheckCLI           = haplocheckCLI,
        mt_dict                 = mt_dict,
        mt_fasta                = mt_fasta,
        mt_fasta_index          = mt_fasta_index,
        mt_amb                  = mt_amb,
        mt_ann                  = mt_ann,
        mt_bwt                  = mt_bwt,
        mt_pac                  = mt_pac,
        mt_sa                   = mt_sa,
        mt_shifted_dict         = mt_shifted_dict,
        mt_shifted_fasta        = mt_shifted_fasta,
        mt_shifted_fasta_index  = mt_shifted_fasta_index,
        mt_shifted_amb          = mt_shifted_amb,
        mt_shifted_ann          = mt_shifted_ann,
        mt_shifted_bwt          = mt_shifted_bwt,
        mt_shifted_pac          = mt_shifted_pac,
        mt_shifted_sa           = mt_shifted_sa,
        compress_output_vcf     = compress_output_vcf,
        shift_back_chain        = shift_back_chain,
        mt_chr_name             = mt_chr_name,
        mt_length               = mt_length,
        mt_nc_start             = mt_nc_start,
        mt_right_pad            = mt_right_pad,
        mt_shift                = mt_shift
    }

    call CoverageAtEveryBase {
      input:
        picard                                           = picard,
        input_bam_regular_ref                            = AlignAndCall.mt_aligned_bam,
        input_bam_regular_ref_index                      = AlignAndCall.mt_aligned_bai,
        input_bam_shifted_ref                            = AlignAndCall.mt_aligned_shifted_bam,
        input_bam_shifted_ref_index                      = AlignAndCall.mt_aligned_shifted_bai,
        shift_back_chain                                 = shift_back_chain,
        control_region_shifted_reference_interval_list   = control_region_shifted_reference_interval_list,
        non_control_region_interval_list                 = non_control_region_interval_list,
        ref_fasta                                        = mt_fasta,
        ref_fasta_index                                  = mt_fasta_index,
        ref_dict                                         = mt_dict,
        shifted_ref_fasta                                = mt_shifted_fasta,
        shifted_ref_fasta_index                          = mt_shifted_fasta_index,
        shifted_ref_dict                                 = mt_shifted_dict,
        mt_length                                        = mt_length,
        shift                                            = mt_shift
    }

    call SplitMultiAllelicSites {
      input:
        gatk            = gatk,
        input_vcf       = AlignAndCall.out_vcf,
        input_vcf_index = AlignAndCall.out_vcf_index,
        base_name       = base_name,
        ref_fasta       = mt_fasta,
        ref_fasta_index = mt_fasta_index,
        ref_dict        = mt_dict
    }
  }

  output {
    Array[File?] contamination_metrics     = AlignAndCall.contamination_metrics
    Array[File] coverage_metrics = AlignAndCall.coverage_metrics
    Array[File] base_level_coverage_metrics  = CoverageAtEveryBase.table
    Array[File] split_vcf                    = SplitMultiAllelicSites.split_vcf
    Array[File] split_vcf_index              = SplitMultiAllelicSites.split_vcf_index
  }
}


# ====================================================================================
# Tasks
# ====================================================================================

task SubsetBamToChrM {
  input {
    String gatk
    File   input_bam
    File   input_bai
    String contig_name
    String basename = basename(basename(input_bam, ".cram"), ".bam")
    File?  ref_fasta
    File?  ref_fasta_index
    File?  ref_dict
  }

  command <<<
    set -euo pipefail

    # 1. 提取所有比对到 chrM 的 Read 及其对应的 Mate (无论 Mate 在哪)
    # -M 确保抓取跨染色体的成对信息
    # 如果是 CRAM，需要提供 -T (ref_fasta)
    samtools view -b \
      ~{if defined(ref_fasta) then "-T ~{ref_fasta}" else ""} \
      -M -o temp_subset.bam \
      ~{input_bam} "~{contig_name}"

    # 2. 此时 temp_subset.bam 可能包含非 chrM 的 Read (即那些 mate 在 chrM 上的 read)
    # 为了保证后续 Picard 工具不报错，我们需要重新排序并建立索引
    # 或者直接进行 RevertSam。但为了稳妥，这里先按 queryname 排序，
    # 这样 RevertSam 处理起来最顺畅。
    
    java -Xmx4G -jar ~{gatk} SortSam \
      -I temp_subset.bam \
      -O ~{basename}.bam \
      --SORT_ORDER queryname
  >>>

  output {
    File output_bam = "~{basename}.bam"
  }
}

task RevertSam {
  input {
    String picard
    File   input_bam
    String basename = basename(input_bam, ".bam")
  }

  meta { description: "Strip alignment while preserving OQ/OA tags and recalibrated qualities." }

  command {
    java -Xmx4G -jar ~{picard} \
      RevertSam \
      INPUT=~{input_bam} \
      OUTPUT_BY_READGROUP=false \
      OUTPUT=~{basename}.bam \
      VALIDATION_STRINGENCY=LENIENT \
      ATTRIBUTE_TO_CLEAR=FT \
      ATTRIBUTE_TO_CLEAR=CO \
      SORT_ORDER=queryname \
      RESTORE_ORIGINAL_QUALITIES=false
  }

  output {
    File unmapped_bam = "~{basename}.bam"
  }
}


task CoverageAtEveryBase {
  input {
    String picard
    File input_bam_regular_ref
    File input_bam_regular_ref_index
    File input_bam_shifted_ref
    File input_bam_shifted_ref_index
    File shift_back_chain
    File control_region_shifted_reference_interval_list
    File non_control_region_interval_list
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File shifted_ref_fasta
    File shifted_ref_fasta_index
    File shifted_ref_dict

    Int mt_length
    Int shift = 8000
  }

  meta { description: "Collect per-base coverage on non-control (unshifted) + control (shifted) regions, then merge." }

  String basename      = basename(input_bam_regular_ref, ".realigned.bam")
  String coverage_file = basename + ".per_base_coverage.tsv"

  command <<<
    set -euo pipefail

    # Collect coverage on the non-control region (using the original reference)
    java -Xmx4G -jar ~{picard} CollectHsMetrics \
      I=~{input_bam_regular_ref} \
      R=~{ref_fasta} \
      PER_BASE_COVERAGE=non_control_region.tsv \
      O=non_control_region.metrics \
      TI=~{non_control_region_interval_list} \
      BI=~{non_control_region_interval_list} \
      COVMAX=20000 \
      SAMPLE_SIZE=1

    # Collect coverage on the control region (using the shifted reference)
    java -Xmx4G -jar ~{picard} CollectHsMetrics \
      I=~{input_bam_shifted_ref} \
      R=~{shifted_ref_fasta} \
      PER_BASE_COVERAGE=control_region_shifted.tsv \
      O=control_region_shifted.metrics \
      TI=~{control_region_shifted_reference_interval_list} \
      BI=~{control_region_shifted_reference_interval_list} \
      COVMAX=20000 \
      SAMPLE_SIZE=1

    # Shift-back control region coordinates and merge with non-control region coverage
    L=~{mt_length}
    S=~{shift}

    # Use awk to map shifted coordinates back and concatenate coverage tables in the correct order
    awk -v L="$L" -v S="$S" -v OFS="\t" '
      BEGIN{ pos_col=0 }
      FNR==1 && FILENAME=="control_region_shifted.tsv" {
        for(i=1;i<=NF;i++){ if($i=="pos"){ pos_col=i } }
        if(pos_col==0){ print "[ERR] pos column not found" > "/dev/stderr"; exit 1 }
        next
      }
      FILENAME=="control_region_shifted.tsv" {
        pos=$pos_col+0;
        pos=((pos-1+S)%L)+1;
        $pos_col=pos;
        if(pos < S){ print > "cr_begin.tsv" }
        else if(pos > S){ print > "cr_end.tsv" }
      }
    ' control_region_shifted.tsv

    head -n 1 non_control_region.tsv > "~{coverage_file}"
    if [[ -s cr_begin.tsv ]]; then cat cr_begin.tsv >> "~{coverage_file}"; fi
    tail -n +2 non_control_region.tsv >> "~{coverage_file}"
    if [[ -s cr_end.tsv ]]; then cat cr_end.tsv >> "~{coverage_file}"; fi
  >>>

  output {
    File table = coverage_file
  }
}


task SplitMultiAllelicSites {
  input {
    String gatk
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File input_vcf
    File input_vcf_index
    String base_name
  }

  String output_vcf       = base_name + ".final.split.vcf"
  String output_vcf_index = output_vcf + ".idx"

  command <<<
    set -e
    java -Xmx4G -jar ~{gatk} LeftAlignAndTrimVariants \
      -R ~{ref_fasta} \
      -V ~{input_vcf} \
      -O ~{output_vcf} \
      --split-multi-allelics \
      --dont-trim-alleles \
      --keep-original-ac
  >>>

  output {
    File split_vcf        = output_vcf
    File split_vcf_index  = output_vcf_index
  }
}