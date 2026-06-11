version 1.0

import "AlignmentPipeline.wdl" as AlignAndMarkDuplicates

workflow AlignAndCall {
  meta { description: "Takes in unmapped bam and outputs VCF of SNP/Indel calls on the mitochondria." }

  input {
    File   unmapped_bam
    Float? autosomal_coverage
    String base_name

    String picard
    String gatk
    String haplocheckCLI

    # original MT reference
    File mt_dict
    File mt_fasta
    File mt_fasta_index
    File mt_amb
    File mt_ann
    File mt_bwt
    File mt_pac
    File mt_sa

    # shifted MT reference (for control region on circular genome)
    File mt_shifted_dict
    File mt_shifted_fasta
    File mt_shifted_fasta_index
    File mt_shifted_amb
    File mt_shifted_ann
    File mt_shifted_bwt
    File mt_shifted_pac
    File mt_shifted_sa

    File shift_back_chain

    # ---- runtime knobs ----
    Boolean compress_output_vcf

    Float?  verifyBamID
    Int?    max_low_het_sites
    String? m2_extra_args
    String? m2_filter_extra_args
    Float?  vaf_filter_threshold
    Float?  f_score_beta
    Int?    max_read_length

    # ==== dynamic region parameters from outer workflow ====
    String? mt_chr_name
    Int?    mt_length
    Int?    mt_nc_start
    Int?    mt_right_pad
    Int?    mt_shift

    # ==== switches ====
    Boolean run_contamination    = false
    Boolean run_blacklist_filter = false

    # optional blacklist
    File?   blacklisted_sites
    File?   blacklisted_sites_index
  }

  # ----------------- Align to original and shifted MT -----------------
  call AlignAndMarkDuplicates.AlignmentPipeline as AlignToMt {
    input:
      input_bam      = unmapped_bam,
      picard         = picard,
      mt_dict        = mt_dict,
      mt_fasta       = mt_fasta,
      mt_fasta_index = mt_fasta_index,
      mt_amb         = mt_amb,
      mt_ann         = mt_ann,
      mt_bwt         = mt_bwt,
      mt_pac         = mt_pac,
      mt_sa          = mt_sa
  }

  call AlignAndMarkDuplicates.AlignmentPipeline as AlignToShiftedMt {
    input:
      input_bam      = unmapped_bam,
      picard         = picard,
      mt_dict        = mt_shifted_dict,
      mt_fasta       = mt_shifted_fasta,
      mt_fasta_index = mt_shifted_fasta_index,
      mt_amb         = mt_shifted_amb,
      mt_ann         = mt_shifted_ann,
      mt_bwt         = mt_shifted_bwt,
      mt_pac         = mt_shifted_pac,
      mt_sa          = mt_shifted_sa
  }

  # coverage metrics
  call CollectWgsMetrics {
    input:
      input_bam       = AlignToMt.mt_aligned_bam,
      input_bam_index = AlignToMt.mt_aligned_bai,
      ref_fasta       = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      picard          = picard,
      read_length     = max_read_length,
      coverage_cap    = 100000
  }

  # ----------------- Call on non-control region (original MT) -----------------
  call M2 as CallMt {
    input:
      input_bam   = AlignToMt.mt_aligned_bam,
      input_bai   = AlignToMt.mt_aligned_bai,
      ref_fasta   = mt_fasta,
      ref_fai     = mt_fasta_index,
      ref_dict    = mt_dict,
      gatk        = gatk,
      compress    = compress_output_vcf,
      m2_extra_args = select_first([m2_extra_args, ""]),
      chr_name    = select_first([mt_chr_name, "chrM"]),
      mt_length   = select_first([mt_length, 16569]),
      nc_start    = select_first([mt_nc_start, 576]),
      right_pad   = select_first([mt_right_pad, 545]),
      shift       = select_first([mt_shift, 8000]),
      call_non_control     = true,
      call_shifted_control = false
  }

  # ----------------- Call on control region (shifted MT) -----------------
  call M2 as CallShiftedMt {
    input:
      input_bam   = AlignToShiftedMt.mt_aligned_bam,
      input_bai   = AlignToShiftedMt.mt_aligned_bai,
      ref_fasta   = mt_shifted_fasta,
      ref_fai     = mt_shifted_fasta_index,
      ref_dict    = mt_shifted_dict,
      gatk        = gatk,
      compress    = compress_output_vcf,
      m2_extra_args = select_first([m2_extra_args, ""]),
      chr_name    = select_first([mt_chr_name, "chrM"]),
      mt_length   = select_first([mt_length, 16569]),
      nc_start    = select_first([mt_nc_start, 576]),
      right_pad   = select_first([mt_right_pad, 545]),
      shift       = select_first([mt_shift, 8000]),
      call_non_control     = false,
      call_shifted_control = true
  }

  # ----------------- Lift over & merge -----------------
  call LiftoverAndCombineVcfs {
    input:
      picard           = picard,
      shifted_vcf      = CallShiftedMt.raw_vcf,
      vcf              = CallMt.raw_vcf,
      ref_fasta        = mt_fasta,
      ref_fasta_index  = mt_fasta_index,
      ref_dict         = mt_dict,
      shift_back_chain = shift_back_chain
      # out_prefix is optional and not required to be passed from outer workflow
  }

  call MergeStats {
    input:
      gatk              = gatk,
      shifted_stats     = CallShiftedMt.stats,
      non_shifted_stats = CallMt.stats
  }

  # ----------------- Initial filtering (always) -----------------
  call Filter as InitialFilter {
    input:
      gatk                      = gatk,
      raw_vcf                   = LiftoverAndCombineVcfs.merged_vcf,
      raw_vcf_index             = LiftoverAndCombineVcfs.merged_vcf_index,
      raw_vcf_stats             = MergeStats.stats,
      base_name                 = base_name,
      ref_fasta                 = mt_fasta,
      ref_fai                   = mt_fasta_index,
      ref_dict                  = mt_dict,
      compress                  = compress_output_vcf,
      m2_extra_filtering_args   = m2_filter_extra_args,
      max_alt_allele_count      = 4,
      vaf_filter_threshold      = 0,
      f_score_beta              = f_score_beta,
      run_blacklist_filter      = run_blacklist_filter,
      blacklisted_sites         = blacklisted_sites,
      blacklisted_sites_index   = blacklisted_sites_index,
      run_contamination         = false
  }

  # ----------------- Optional contamination path -----------------
  if (run_contamination) {
    call SplitMultiAllelicsAndRemoveNonPassSites {
      input:
        gatk             = gatk,
        ref_fasta        = mt_fasta,
        ref_fai          = mt_fasta_index,
        ref_dict         = mt_dict,
        filtered_vcf     = InitialFilter.filtered_vcf,
        filtered_vcf_idx = InitialFilter.filtered_vcf_idx
    }

    call GetContamination {
      input:
        haplocheckCLI = haplocheckCLI,
        input_vcf     = SplitMultiAllelicsAndRemoveNonPassSites.vcf_for_haplochecker
    }

    call Filter as FilterContamination {
      input:
        gatk                    = gatk,
        raw_vcf                 = InitialFilter.filtered_vcf,
        raw_vcf_index           = InitialFilter.filtered_vcf_idx,
        raw_vcf_stats           = MergeStats.stats,
        run_contamination       = true,
        hasContamination        = GetContamination.hasContamination,
        contamination_major     = GetContamination.major_level,
        contamination_minor     = GetContamination.minor_level,
        verifyBamID             = verifyBamID,
        base_name               = base_name,
        ref_fasta               = mt_fasta,
        ref_fai                 = mt_fasta_index,
        ref_dict                = mt_dict,
        compress                = compress_output_vcf,
        m2_extra_filtering_args = m2_filter_extra_args,
        max_alt_allele_count    = 4,
        vaf_filter_threshold    = vaf_filter_threshold,
        f_score_beta            = f_score_beta,
        run_blacklist_filter    = run_blacklist_filter,
        blacklisted_sites       = blacklisted_sites,
        blacklisted_sites_index = blacklisted_sites_index
    }
  }

  # choose filtered VCF
  File  chosen_vcf      = select_first([FilterContamination.filtered_vcf,     InitialFilter.filtered_vcf])
  File  chosen_vcf_idx  = select_first([FilterContamination.filtered_vcf_idx, InitialFilter.filtered_vcf_idx])

  if ( defined(autosomal_coverage) ) {
    call FilterNuMTs {
      input:
        gatk               = gatk,
        filtered_vcf       = chosen_vcf,
        filtered_vcf_index = chosen_vcf_idx,
        ref_fasta          = mt_fasta,
        ref_fai            = mt_fasta_index,
        ref_dict           = mt_dict,
        autosomal_coverage = autosomal_coverage,
        compress           = compress_output_vcf
    }
  }

  File low_het_vcf       = select_first([FilterNuMTs.numt_filtered_vcf,     chosen_vcf])
  File low_het_vcf_index = select_first([FilterNuMTs.numt_filtered_vcf_idx, chosen_vcf_idx])

  call FilterLowHetSites {
    input:
      gatk               = gatk,
      filtered_vcf       = low_het_vcf,
      filtered_vcf_index = low_het_vcf_index,
      ref_fasta          = mt_fasta,
      ref_fai            = mt_fasta_index,
      ref_dict           = mt_dict,
      max_low_het_sites  = max_low_het_sites,
      compress           = compress_output_vcf,
      base_name          = base_name
  }

  output {
    File  mt_aligned_bam         = AlignToMt.mt_aligned_bam
    File  mt_aligned_bai         = AlignToMt.mt_aligned_bai
    File  mt_aligned_shifted_bam = AlignToShiftedMt.mt_aligned_bam
    File  mt_aligned_shifted_bai = AlignToShiftedMt.mt_aligned_bai
    File  duplicate_metrics      = AlignToMt.duplicate_metrics
    File  coverage_metrics       = CollectWgsMetrics.metrics
    File  theoretical_sensitivity_metrics = CollectWgsMetrics.theoretical_sensitivity
    Int   mean_coverage          = CollectWgsMetrics.mean_coverage

    File?  contamination_metrics       = GetContamination.contamination_file
    String? major_haplogroup           = GetContamination.major_hg
    Float   contamination              = select_first([FilterContamination.contamination, 0.0])
    File?   input_vcf_for_haplochecker = SplitMultiAllelicsAndRemoveNonPassSites.vcf_for_haplochecker

    File out_vcf       = FilterLowHetSites.final_filtered_vcf
    File out_vcf_index = FilterLowHetSites.final_filtered_vcf_idx
  }
}


workflow AlignAndCallPrealigned {
  meta {
    description: "mtSwirl-like Round 2 caller. Takes prealigned, duplicate-marked, chrM-assigned BAMs from competitive chrM+NUMT self-reference alignment and outputs final filtered VCF."
  }

  input {
    File   standard_bam
    File   standard_bai
    File   shifted_bam
    File   shifted_bai
    Float? autosomal_coverage
    String base_name

    String picard
    String gatk
    String haplocheckCLI

    # Standard self-reference. For full mtSwirl-like behavior this can contain chrM + NUMT contigs.
    File mt_dict
    File mt_fasta
    File mt_fasta_index

    # Shifted self-reference. For full mtSwirl-like behavior this can contain shifted chrM + NUMT contigs.
    File mt_shifted_dict
    File mt_shifted_fasta
    File mt_shifted_fasta_index

    File shift_back_chain

    Boolean compress_output_vcf

    Float?  verifyBamID
    Int?    max_low_het_sites
    String? m2_extra_args
    String? m2_filter_extra_args
    Float?  vaf_filter_threshold
    Float?  f_score_beta
    Int?    max_read_length

    String? mt_chr_name
    Int?    mt_length
    Int?    mt_nc_start
    Int?    mt_right_pad
    Int?    mt_shift

    Boolean run_contamination    = false
    Boolean run_blacklist_filter = false

    File?   blacklisted_sites
    File?   blacklisted_sites_index

    # Optional force-call sites in the same coordinates/dictionary as mt_fasta.
    File?   force_call_sites_vcf
    File?   force_call_sites_vcf_index
  }

  call CollectWgsMetrics {
    input:
      input_bam       = standard_bam,
      input_bam_index = standard_bai,
      ref_fasta       = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      picard          = picard,
      read_length     = max_read_length,
      coverage_cap    = 100000
  }

  call M2 as CallMt {
    input:
      input_bam   = standard_bam,
      input_bai   = standard_bai,
      ref_fasta   = mt_fasta,
      ref_fai     = mt_fasta_index,
      ref_dict    = mt_dict,
      gatk        = gatk,
      compress    = compress_output_vcf,
      m2_extra_args = select_first([m2_extra_args, ""]),
      chr_name    = select_first([mt_chr_name, "chrM"]),
      mt_length   = select_first([mt_length, 16569]),
      nc_start    = select_first([mt_nc_start, 576]),
      right_pad   = select_first([mt_right_pad, 545]),
      shift       = select_first([mt_shift, 8000]),
      call_non_control     = true,
      call_shifted_control = false,
      force_call_sites_vcf = force_call_sites_vcf,
      force_call_sites_vcf_index = force_call_sites_vcf_index,
      use_mate_contig_read_filters = false
  }

  call M2 as CallShiftedMt {
    input:
      input_bam   = shifted_bam,
      input_bai   = shifted_bai,
      ref_fasta   = mt_shifted_fasta,
      ref_fai     = mt_shifted_fasta_index,
      ref_dict    = mt_shifted_dict,
      gatk        = gatk,
      compress    = compress_output_vcf,
      m2_extra_args = select_first([m2_extra_args, ""]),
      chr_name    = select_first([mt_chr_name, "chrM"]),
      mt_length   = select_first([mt_length, 16569]),
      nc_start    = select_first([mt_nc_start, 576]),
      right_pad   = select_first([mt_right_pad, 545]),
      shift       = select_first([mt_shift, 8000]),
      call_non_control     = false,
      call_shifted_control = true,
      force_call_sites_vcf = force_call_sites_vcf,
      force_call_sites_vcf_index = force_call_sites_vcf_index,
      use_mate_contig_read_filters = false
  }

  call LiftoverAndCombineVcfs {
    input:
      picard           = picard,
      shifted_vcf      = CallShiftedMt.raw_vcf,
      vcf              = CallMt.raw_vcf,
      ref_fasta        = mt_fasta,
      ref_fasta_index  = mt_fasta_index,
      ref_dict         = mt_dict,
      shift_back_chain = shift_back_chain
  }

  call MergeStats {
    input:
      gatk              = gatk,
      shifted_stats     = CallShiftedMt.stats,
      non_shifted_stats = CallMt.stats
  }

  call Filter as InitialFilter {
    input:
      gatk                      = gatk,
      raw_vcf                   = LiftoverAndCombineVcfs.merged_vcf,
      raw_vcf_index             = LiftoverAndCombineVcfs.merged_vcf_index,
      raw_vcf_stats             = MergeStats.stats,
      base_name                 = base_name,
      ref_fasta                 = mt_fasta,
      ref_fai                   = mt_fasta_index,
      ref_dict                  = mt_dict,
      compress                  = compress_output_vcf,
      m2_extra_filtering_args   = m2_filter_extra_args,
      max_alt_allele_count      = 4,
      vaf_filter_threshold      = 0,
      f_score_beta              = f_score_beta,
      run_blacklist_filter      = run_blacklist_filter,
      blacklisted_sites         = blacklisted_sites,
      blacklisted_sites_index   = blacklisted_sites_index,
      run_contamination         = false
  }

  if (run_contamination) {
    call SplitMultiAllelicsAndRemoveNonPassSites {
      input:
        gatk             = gatk,
        ref_fasta        = mt_fasta,
        ref_fai          = mt_fasta_index,
        ref_dict         = mt_dict,
        filtered_vcf     = InitialFilter.filtered_vcf,
        filtered_vcf_idx = InitialFilter.filtered_vcf_idx
    }

    call GetContamination {
      input:
        haplocheckCLI = haplocheckCLI,
        input_vcf     = SplitMultiAllelicsAndRemoveNonPassSites.vcf_for_haplochecker
    }

    call Filter as FilterContamination {
      input:
        gatk                    = gatk,
        raw_vcf                 = InitialFilter.filtered_vcf,
        raw_vcf_index           = InitialFilter.filtered_vcf_idx,
        raw_vcf_stats           = MergeStats.stats,
        run_contamination       = true,
        hasContamination        = GetContamination.hasContamination,
        contamination_major     = GetContamination.major_level,
        contamination_minor     = GetContamination.minor_level,
        verifyBamID             = verifyBamID,
        base_name               = base_name,
        ref_fasta               = mt_fasta,
        ref_fai                 = mt_fasta_index,
        ref_dict                = mt_dict,
        compress                = compress_output_vcf,
        m2_extra_filtering_args = m2_filter_extra_args,
        max_alt_allele_count    = 4,
        vaf_filter_threshold    = vaf_filter_threshold,
        f_score_beta            = f_score_beta,
        run_blacklist_filter    = run_blacklist_filter,
        blacklisted_sites       = blacklisted_sites,
        blacklisted_sites_index = blacklisted_sites_index
    }
  }

  File  chosen_vcf      = select_first([FilterContamination.filtered_vcf,     InitialFilter.filtered_vcf])
  File  chosen_vcf_idx  = select_first([FilterContamination.filtered_vcf_idx, InitialFilter.filtered_vcf_idx])

  if ( defined(autosomal_coverage) ) {
    call FilterNuMTs {
      input:
        gatk               = gatk,
        filtered_vcf       = chosen_vcf,
        filtered_vcf_index = chosen_vcf_idx,
        ref_fasta          = mt_fasta,
        ref_fai            = mt_fasta_index,
        ref_dict           = mt_dict,
        autosomal_coverage = autosomal_coverage,
        compress           = compress_output_vcf
    }
  }

  File low_het_vcf       = select_first([FilterNuMTs.numt_filtered_vcf,     chosen_vcf])
  File low_het_vcf_index = select_first([FilterNuMTs.numt_filtered_vcf_idx, chosen_vcf_idx])

  call FilterLowHetSites {
    input:
      gatk               = gatk,
      filtered_vcf       = low_het_vcf,
      filtered_vcf_index = low_het_vcf_index,
      ref_fasta          = mt_fasta,
      ref_fai            = mt_fasta_index,
      ref_dict           = mt_dict,
      max_low_het_sites  = max_low_het_sites,
      compress           = compress_output_vcf,
      base_name          = base_name
  }

  output {
    File  mt_aligned_bam         = standard_bam
    File  mt_aligned_bai         = standard_bai
    File  mt_aligned_shifted_bam = shifted_bam
    File  mt_aligned_shifted_bai = shifted_bai
    File  duplicate_metrics      = standard_bam
    File  coverage_metrics       = CollectWgsMetrics.metrics
    File  theoretical_sensitivity_metrics = CollectWgsMetrics.theoretical_sensitivity
    Int   mean_coverage          = CollectWgsMetrics.mean_coverage

    File?  contamination_metrics       = GetContamination.contamination_file
    String? major_haplogroup           = GetContamination.major_hg
    Float   contamination              = select_first([FilterContamination.contamination, 0.0])
    File?   input_vcf_for_haplochecker = SplitMultiAllelicsAndRemoveNonPassSites.vcf_for_haplochecker

    File out_vcf       = FilterLowHetSites.final_filtered_vcf
    File out_vcf_index = FilterLowHetSites.final_filtered_vcf_idx
  }
}

# ===================== TASKS =====================

task CollectWgsMetrics {
  input {
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    String picard
    Int? read_length
    Int? coverage_cap
  }

  Int read_length_for_optimization = select_first([read_length, 151])

  command <<<
  set -e

  java -Xms2000m -jar ~{picard} \
    CollectWgsMetrics \
    INPUT=~{input_bam} \
    VALIDATION_STRINGENCY=SILENT \
    REFERENCE_SEQUENCE=~{ref_fasta} \
    OUTPUT=metrics.txt \
    USE_FAST_ALGORITHM=true \
    READ_LENGTH=~{read_length_for_optimization} \
    ~{"COVERAGE_CAP=" + coverage_cap} \
    INCLUDE_BQ_HISTOGRAM=true \
    THEORETICAL_SENSITIVITY_OUTPUT=theoretical_sensitivity.txt

  # Parse MEAN_COVERAGE without R:
  awk -F'\t' '
    BEGIN{col=-1}
    $1=="GENOME_TERRITORY"{
      for(i=1;i<=NF;i++) if($i=="MEAN_COVERAGE") col=i; next
    }
    col>0 && $1!~/^#/ && $1!~/^GENOME_TERRITORY/{
      val=$col+0; printf("%d\n", val); exit
    }
  ' metrics.txt > mean_coverage.txt
  >>>


  output {
    File metrics                 = "metrics.txt"
    File theoretical_sensitivity = "theoretical_sensitivity.txt"
    Int  mean_coverage           = read_int("mean_coverage.txt")
  }
}

task M2 {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File input_bam
    File input_bai
    String gatk
    Int? max_reads_per_alignment_start
    String? m2_extra_args
    Boolean? make_bamout
    Boolean compress

    # dynamic region
    String chr_name = "chrM"
    Int    mt_length = 16569
    Int    nc_start  = 576
    Int    right_pad = 545
    Int    shift     = 8000
    Boolean call_non_control     = true
    Boolean call_shifted_control = false

    # mtSwirl-like Round 2 additions.
    # Provide a VCF of consensus-basis homoplasmies in self-reference coordinates if available.
    File? force_call_sites_vcf
    File? force_call_sites_vcf_index

    # Existing GATK mitochondrial workflow used mate-contig filters.
    # For mtSwirl-like preassigned chrM reads, set this to false because reads whose mate
    # was assigned to NUMT should not be discarded after competitive realignment.
    Boolean use_mate_contig_read_filters = true
  }

  Int    max_reads_per_alignment_start_arg = select_first([max_reads_per_alignment_start, 75])
  String output_vcf        = "raw" + (if compress then ".vcf.gz" else ".vcf")
  String output_vcf_index  = output_vcf + (if compress then ".tbi" else ".idx")

  command <<<
    set -euo pipefail

    touch bamout.bam

    CHR="~{chr_name}"
    L=~{mt_length}
    NC_START=~{nc_start}
    RIGHT_PAD=~{right_pad}
    NC_END=$(( L - RIGHT_PAD ))

    if (( NC_END < NC_START )); then
      echo "[M2] ERROR: NC_END(${NC_END}) < NC_START(${NC_START}); L=${L}, RIGHT_PAD=${RIGHT_PAD}" >&2
      exit 1
    fi

    REGION=""
    if [[ "~{call_non_control}" == "true" ]]; then
      REGION="${CHR}:${NC_START}-${NC_END}"
    fi
    if [[ "~{call_shifted_control}" == "true" ]]; then
      S=~{shift}
      start=$(( ((NC_END   - S) % L + L) % L + 1 ))
      end=$((   ((NC_START - 1 - S) % L + L) % L + 1 ))
      REGION="${CHR}:${start}-${end}"
    fi

    echo "[M2] CHR=$CHR L=$L NC_START=$NC_START NC_END=$NC_END REGION=$REGION" >&2

    # Localize force-call index if provided.
    ~{if defined(force_call_sites_vcf_index) then "ls " + force_call_sites_vcf_index + " >/dev/null" else ":"}

    java -Xmx4G -jar ~{gatk} Mutect2 \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      ~{if use_mate_contig_read_filters then "--read-filter MateOnSameContigOrNoMappedMateReadFilter --read-filter MateUnmappedAndUnmappedReadFilter" else ""} \
      ~{if defined(force_call_sites_vcf) then "--alleles " + force_call_sites_vcf else ""} \
      -O ~{output_vcf} \
      ~{true='--bam-output bamout.bam' false='' make_bamout} \
      ~{m2_extra_args} \
      --annotation StrandBiasBySample \
      --mitochondria-mode \
      --max-reads-per-alignment-start ~{max_reads_per_alignment_start_arg} \
      --max-mnp-distance 0 \
      -L "${REGION}"
  >>>

  output {
    File raw_vcf       = output_vcf
    File raw_vcf_idx   = output_vcf_index
    File stats         = output_vcf + ".stats"
    File output_bamOut = "bamout.bam"
  }
}

task LiftoverAndCombineVcfs {
  input {
    String picard
    File shifted_vcf
    File vcf
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File shift_back_chain
    String? out_prefix
  }

  String base = select_first([out_prefix, basename(shifted_vcf, ".vcf")])

  command <<<
    set -e
    java -Xms2000m -jar ~{picard} LiftoverVcf \
      I=~{shifted_vcf} \
      O=~{base}.shifted_back.vcf \
      R=~{ref_fasta} \
      CHAIN=~{shift_back_chain} \
      REJECT=~{base}.rejected.vcf

    java -Xms2000m -jar ~{picard} MergeVcfs \
      I=~{base}.shifted_back.vcf \
      I=~{vcf} \
      O=~{base}.merged.vcf
  >>>

  output {
    File rejected_vcf      = base + ".rejected.vcf"
    File merged_vcf        = base + ".merged.vcf"
    File merged_vcf_index  = base + ".merged.vcf.idx"
  }
}

task MergeStats {
  input {
    String gatk
    File shifted_stats
    File non_shifted_stats
  }

  command <<<
    set -e
    java -Xmx4G -jar ~{gatk} MergeMutectStats \
      --stats ~{shifted_stats} \
      --stats ~{non_shifted_stats} \
      -O raw.combined.stats
  >>>

  output { File stats = "raw.combined.stats" }
}

task Filter {
  input {
    String gatk
    File ref_fasta
    File ref_fai
    File ref_dict
    File raw_vcf
    File raw_vcf_index
    File raw_vcf_stats
    Boolean compress
    Float? vaf_cutoff
    String base_name

    String? m2_extra_filtering_args
    Int max_alt_allele_count
    Float? autosomal_coverage
    Float? vaf_filter_threshold
    Float? f_score_beta

    Boolean run_contamination
    String? hasContamination
    Float? contamination_major
    Float? contamination_minor
    Float? verifyBamID

    Boolean run_blacklist_filter
    File?   blacklisted_sites
    File?   blacklisted_sites_index
  }

  String output_vcf       = base_name + (if compress then ".vcf.gz" else ".vcf")
  String output_vcf_index = output_vcf + (if compress then ".tbi" else ".idx")
  String pre_vcf          = "pre_blacklist" + (if compress then ".vcf.gz" else ".vcf")
  String pre_vcf_index    = pre_vcf + (if compress then ".tbi" else ".idx")

  Float hc_contamination  = if run_contamination && hasContamination == "YES" then (if contamination_major == 0.0 then contamination_minor else 1.0 - contamination_major) else 0.0
  Float max_contamination = if defined(verifyBamID) && verifyBamID > hc_contamination then verifyBamID else hc_contamination

  command <<<
    set -e

    java -Xmx4G -jar ~{gatk} FilterMutectCalls \
      -V ~{raw_vcf} \
      -R ~{ref_fasta} \
      -O ~{pre_vcf} \
      --stats ~{raw_vcf_stats} \
      ~{m2_extra_filtering_args} \
      --max-alt-allele-count ~{max_alt_allele_count} \
      --mitochondria-mode \
      ~{"--min-allele-fraction " + vaf_filter_threshold} \
      ~{"--f-score-beta " + f_score_beta} \
      ~{"--contamination-estimate " + max_contamination}

    if [[ "~{run_blacklist_filter}" == "true" && -s "~{blacklisted_sites}" ]]; then
      java -Xmx4G -jar ~{gatk} VariantFiltration \
        -V ~{pre_vcf} \
        -O ~{output_vcf} \
        --apply-allele-specific-filters \
        --mask ~{blacklisted_sites} \
        --mask-name "blacklisted_site"
    else
      mv ~{pre_vcf} ~{output_vcf}
      if [[ -f "~{pre_vcf_index}" ]]; then
        mv ~{pre_vcf_index} ~{output_vcf_index}
      fi
    fi
  >>>

  output {
    File  filtered_vcf     = output_vcf
    File  filtered_vcf_idx = output_vcf_index
    Float contamination    = hc_contamination
  }
}

task SplitMultiAllelicsAndRemoveNonPassSites {
  input {
    String gatk
    File ref_fasta
    File ref_fai
    File ref_dict
    File filtered_vcf
    File filtered_vcf_idx
  }

  String basename   = basename(filtered_vcf, ".vcf.gz")
  String output_vcf = basename + ".splitAndPassOnly.vcf"

  command <<<
    set -e
    java -Xmx4G -jar ~{gatk} LeftAlignAndTrimVariants \
      -R ~{ref_fasta} \
      -V ~{filtered_vcf} \
      -O split.vcf \
      --split-multi-allelics \
      --dont-trim-alleles \
      --keep-original-ac

    java -Xmx4G -jar ~{gatk} SelectVariants \
      -V split.vcf \
      -O ~{output_vcf} \
      --exclude-filtered
  >>>

  output { File vcf_for_haplochecker = output_vcf }
}

task GetContamination {
  input {
    String haplocheckCLI
    File input_vcf
  }

  String basename    = basename(input_vcf, ".splitAndPassOnly.vcf")
  String output_file = basename + ".haplocheck_contamination.txt"

  command <<<
    set -e
    PARENT_DIR="$(dirname "~{input_vcf}")"
    java -Xmx4G -jar ~{haplocheckCLI} "${PARENT_DIR}"

    sed 's/\"//g' output > ~{output_file}
    grep "SampleID" ~{output_file} > headers

    FORMAT_ERROR="Bad contamination file format"
    [[ $(awk '{print $2}'  headers) == "Contamination"     ]] || { echo $FORMAT_ERROR; exit 1; }
    [[ $(awk '{print $6}'  headers) == "HgMajor"           ]] || { echo $FORMAT_ERROR; exit 1; }
    [[ $(awk '{print $8}'  headers) == "HgMinor"           ]] || { echo $FORMAT_ERROR; exit 1; }
    [[ $(awk '{print $14}' headers) == "MeanHetLevelMajor" ]] || { echo $FORMAT_ERROR; exit 1; }
    [[ $(awk '{print $15}' headers) == "MeanHetLevelMinor" ]] || { echo $FORMAT_ERROR; exit 1; }

    grep -v "SampleID" ~{output_file} > output-data
    awk -F "\t" '{print $2}'  output-data > contamination.txt
    awk -F "\t" '{print $6}'  output-data > major_hg.txt
    awk -F "\t" '{print $8}'  output-data > minor_hg.txt
    awk -F "\t" '{print $14}' output-data > mean_het_major.txt
    awk -F "\t" '{print $15}' output-data > mean_het_minor.txt
  >>>

  output {
    File  contamination_file = output_file
    String hasContamination  = read_string("contamination.txt")
    String major_hg          = read_string("major_hg.txt")
    String minor_hg          = read_string("minor_hg.txt")
    Float  major_level       = read_float("mean_het_major.txt")
    Float  minor_level       = read_float("mean_het_minor.txt")
  }
}

task FilterNuMTs {
  input {
    String gatk
    File ref_fasta
    File ref_fai
    File ref_dict
    File filtered_vcf
    File filtered_vcf_index
    Float? autosomal_coverage
    Boolean compress
  }

  String basename   = basename(filtered_vcf, ".vcf")
  String output_vcf = basename + ".numt" + (if compress then ".vcf.gz" else ".vcf")
  String output_vcf_index = output_vcf + (if compress then ".tbi" else ".idx")

  command <<<
    set -e
    java -Xmx4G -jar ~{gatk} NuMTFilterTool \
      -R ~{ref_fasta} \
      -V ~{filtered_vcf} \
      -O ~{output_vcf} \
      --autosomal-coverage ~{autosomal_coverage}
  >>>

  output {
    File numt_filtered_vcf     = output_vcf
    File numt_filtered_vcf_idx = output_vcf_index
  }
}

task FilterLowHetSites {
  input {
    String gatk
    File ref_fasta
    File ref_fai
    File ref_dict
    File filtered_vcf
    File filtered_vcf_index
    String base_name
    Int? max_low_het_sites
    Boolean compress
  }

  String output_vcf       = base_name + ".final" + (if compress then ".vcf.gz" else ".vcf")
  String output_vcf_index = output_vcf + (if compress then ".tbi" else ".idx")
  Int    max_sites        = select_first([max_low_het_sites, 1000])

  command <<<
    set -e
    java -Xmx4G -jar ~{gatk} MTLowHeteroplasmyFilterTool \
      -R ~{ref_fasta} \
      -V ~{filtered_vcf} \
      -O ~{output_vcf} \
      --max-allowed-low-hets ~{max_sites}
  >>>

  output {
    File final_filtered_vcf     = output_vcf
    File final_filtered_vcf_idx = output_vcf_index
  }
}
