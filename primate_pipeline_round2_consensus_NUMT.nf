#!/usr/bin/env nextflow
nextflow.enable.dsl=2

if (!params.sample_tsv) {
    error "Missing required input: --sample_tsv sample.tsv"
}
if (!params.round1_outdir) {
    error "Missing required input: --round1_outdir /path/to/round1_output_dir"
}
if (!params.ref_dir) {
    error "Missing required input: --ref_dir /path/to/mt_reference_dir"
}
if (!params.wdl_script_round2) {
    error "Missing required input: --wdl_script_round2 /path/to/MitochondriaPipeline.wdl"
}
if (!params.cromwell_jar) {
    error "Missing required input: --cromwell_jar /path/to/cromwell.jar"
}
if (!params.picard_jar) {
    error "Missing required input: --picard_jar /path/to/picard.jar"
}
if (!params.gatk_jar) {
    error "Missing required input: --gatk_jar /path/to/gatk.jar"
}
if (!params.haplocheck_jar) {
    error "Missing required input: --haplocheck_jar /path/to/haplocheck.jar"
}

// Optional: final round2 VCF name produced by the WDL inside Cromwell outputs.
// The lift-back step searches recursively under cromwell-executions for this file.
params.round2_final_vcf_suffix = params.round2_final_vcf_suffix ?: '.round2.consensus.final.split.vcf'

// Round1 original NUMT FASTA published by PREPARE_DECOY_REFERENCE.
// Expected path by default:
//   <round1_outdir>/<sample_id>/<round1_numt_subdir>/<sample_id><round1_numt_suffix>
params.round1_numt_subdir = params.round1_numt_subdir ?: 'numt_decoy_ref'
params.round1_numt_suffix = params.round1_numt_suffix ?: '.original_numt.fa'
// Round1 decoy-coordinate NUMT VCF published by CALL_NUMT_VARIANTS_DECOY.
// Expected path by default:
//   <round1_outdir>/<sample_id>/round_1/<round1_numt_vcf_subdir>/<sample_id><round1_nuc_vcf_suffix>
params.round1_numt_vcf_subdir = params.round1_numt_vcf_subdir ?: 'numt_decoy_variant_calling'
params.round1_nuc_vcf_suffix = params.round1_nuc_vcf_suffix ?: '.numt_decoy.raw.vcf.gz'
params.strict_numt_ref = params.strict_numt_ref ?: true

log.info """
ROUND2 CONSENSUS mtDNA PIPELINE
===============================
Sample TSV:             ${params.sample_tsv}
Round1 output dir:      ${params.round1_outdir}
Round1 BAM subdir:      ${params.round1_bam_subdir}
Round1 VCF subdir:      ${params.round1_vcf_subdir}
Output dir:             ${params.outdir}
Consensus filter expr:  ${params.consensus_filter_expr}
mt ref dir:             ${params.ref_dir}
Round1 NUMT subdir:     ${params.round1_numt_subdir}
Round1 NUMT suffix:     ${params.round1_numt_suffix}
Round1 NUMT VCF subdir: ${params.round1_numt_vcf_subdir}
Round1 NUMT VCF suffix: ${params.round1_nuc_vcf_suffix}
Strict NUMT ref:        ${params.strict_numt_ref}
mt shift:               ${params.mt_shift}
interval padding:       ${params.mt_interval_padding}
WDL ref files:          generated inside this NF
WDL script:             ${params.wdl_script_round2}
===============================
"""

ch_samples = Channel.fromPath(params.sample_tsv)
    .splitCsv(header: false, sep: '\t')
    .filter { row -> row.size() >= 3 && row[0]?.trim() && !row[0].startsWith('#') }
    .map { row ->
        def sample_id = row[0].trim()
        def species_name = row[1].trim()
        def ref_name = row[2].trim()

        def round2_dir = file("${params.outdir}/${sample_id}/round_2")
        def original_coords_dir = file("${params.outdir}/${sample_id}/round_2_variant_calling_original_coords")

        if (round2_dir.exists() && original_coords_dir.exists()) {
            log.info "SKIP completed sample ${sample_id}: ${round2_dir} and ${original_coords_dir} already exist"
            return null
        }

        def meta = [id: sample_id]
        tuple(meta, species_name, ref_name)
    }
    .filter { it != null }

ch_cromwell_conf = file("${baseDir}/cromwell.conf")

workflow {
    FIND_ROUND1_OUTPUTS(ch_samples)

    ch_round1_files = FIND_ROUND1_OUTPUTS.out.round1_files
    ch_vcf_inputs = ch_round1_files.map { meta, species_name, ref_name, round1_vcf, round1_vcf_tbi, round1_bam, original_numt_fa, original_numt_fai, round1_nuc_vcf, round1_nuc_vcf_tbi ->
        tuple(meta, species_name, ref_name, round1_vcf, original_numt_fa, original_numt_fai, round1_nuc_vcf, round1_nuc_vcf_tbi)
    }
    ch_bam_inputs = ch_round1_files.map { meta, species_name, ref_name, round1_vcf, round1_vcf_tbi, round1_bam, original_numt_fa, original_numt_fai, round1_nuc_vcf, round1_nuc_vcf_tbi ->
        tuple(meta, species_name, ref_name, round1_bam)
    }

    BUILD_CONSENSUS_REFERENCE(ch_vcf_inputs)
    ch_bam_consensus = ch_bam_inputs.join(BUILD_CONSENSUS_REFERENCE.out.consensus_ref, by: [0,1,2])
    REALIGN_TO_CONSENSUS_ASSIGNED_BAMS(ch_bam_consensus)
    GENERATE_BAM_TSV(REALIGN_TO_CONSENSUS_ASSIGNED_BAMS.out.assigned_bams)
    ch_json_inputs = GENERATE_BAM_TSV.out.tsv.join(BUILD_CONSENSUS_REFERENCE.out.consensus_ref, by: [0,1,2])
    GENERATE_WDL_JSON_ROUND2(ch_json_inputs)
    RUN_WDL_VARIANT_CALLING(GENERATE_WDL_JSON_ROUND2.out.json, ch_cromwell_conf)

    // Convert round2 VCF calls from sample-consensus coordinates back to the
    // original mt reference coordinates, using the exact SNV/indel set that
    // was used to build the consensus reference.
    ch_liftback_inputs = RUN_WDL_VARIANT_CALLING.out.cromwell_out.join(BUILD_CONSENSUS_REFERENCE.out.consensus_ref, by: [0,1,2])
    LIFTBACK_ROUND2_VCF_TO_ORIGINAL(ch_liftback_inputs)
}

workflow.onComplete {
    if (workflow.success) {
        log.info "Round2 consensus pipeline completed successfully. Output files are in: ${params.outdir}"
    } else {
        log.error "Round2 consensus pipeline finished with errors. Check .nextflow.log and failed task work dirs."
    }
}

process FIND_ROUND1_OUTPUTS {
    tag { "Find round1 outputs for ${meta.id}" }
    label 'generation_related'
    publishDir "${params.outdir}/${meta.id}/round_2/round2_inputs_from_round1", mode: 'copy'

    input:
    tuple val(meta), val(species_name), val(ref_name)

    output:
    tuple val(meta), val(species_name), val(ref_name),
          path("${meta.id}.round1.input.vcf.gz"),
          path("${meta.id}.round1.input.vcf.gz.tbi"),
          path("${meta.id}.round1.input.bam"),
          path("${meta.id}.round1.original_numt.fa"),
          path("${meta.id}.round1.original_numt.fa.fai"),
          path("${meta.id}.round1.nuc.input.vcf.gz"),
          path("${meta.id}.round1.nuc.input.vcf.gz.tbi"),
          emit: round1_files

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    ROUND1_ROOT="${params.round1_outdir}"
    SAMPLE_ID="${meta.id}"
    SAMPLE_DIR="\${ROUND1_ROOT}/\${SAMPLE_ID}"

    BAM_ROOT="\${SAMPLE_DIR}/round_1/${params.round1_bam_subdir}"
    VCF_ROOT="\${SAMPLE_DIR}/${params.round1_vcf_subdir}"
    NUMT_ROOT="\${SAMPLE_DIR}/round_1/${params.round1_numt_subdir}"
    NUMT_VCF_ROOT="\${SAMPLE_DIR}/round_1/${params.round1_numt_vcf_subdir}"

    [[ -d "\${SAMPLE_DIR}" ]] || {
        echo "ERROR: Missing round1 sample directory: \${SAMPLE_DIR}" >&2
        exit 1
    }

    [[ -d "\${BAM_ROOT}" ]] || {
        echo "ERROR: Missing round1 BAM directory: \${BAM_ROOT}" >&2
        exit 1
    }

    [[ -d "\${VCF_ROOT}" ]] || {
        echo "ERROR: Missing round1 VCF root directory: \${VCF_ROOT}" >&2
        exit 1
    }

    [[ -d "\${NUMT_ROOT}" ]] || {
        echo "ERROR: Missing round1 NUMT FASTA directory: \${NUMT_ROOT}" >&2
        exit 1
    }

    [[ -d "\${NUMT_VCF_ROOT}" ]] || {
        echo "[WARN] Missing round1 decoy NUMT VCF directory: \${NUMT_VCF_ROOT}; consensus NUMT will equal original NUMT if no fallback VCF is found" >&2
    }

    ########################################
    # 1) Use candidate reads BAM.
    #    No BAM index is required for samtools collate / fastq.
    ########################################
    bam="\${BAM_ROOT}/\${SAMPLE_ID}${params.round1_bam_suffix}"

    [[ -s "\${bam}" ]] || {
        echo "ERROR: Missing round1 candidate BAM: \${bam}" >&2
        exit 1
    }

    ln -sf "\${bam}" "\${SAMPLE_ID}.round1.input.bam"

    ########################################
    # 2) Use sample-specific original NUMT FASTA from round1 PREPARE_DECOY_REFERENCE.
    #    Round 2 will convert this to consensus NUMT using the Round 1 NUMT VCF.
    ########################################
    numt_fa="\${NUMT_ROOT}/\${SAMPLE_ID}${params.round1_numt_suffix}"
    numt_fai="\${numt_fa}.fai"

    [[ -e "\${numt_fa}" ]] || {
        echo "ERROR: Missing round1 original NUMT FASTA: \${numt_fa}" >&2
        exit 1
    }

    if [[ ! -s "\${numt_fa}" ]]; then
        echo "[INFO] Round1 original NUMT FASTA is empty for \${SAMPLE_ID}; continuing with chrM-only self-reference" >&2
        : > "\${SAMPLE_ID}.round1.original_numt.fa"
        : > "\${SAMPLE_ID}.round1.original_numt.fa.fai"
    elif [[ ! -s "\${numt_fai}" ]]; then
        echo "[WARN] Missing original NUMT FASTA index; creating a local index for \${numt_fa}" >&2
        cp "\${numt_fa}" "\${SAMPLE_ID}.round1.original_numt.fa"
        samtools faidx "\${SAMPLE_ID}.round1.original_numt.fa"
    else
        ln -sf "\${numt_fa}" "\${SAMPLE_ID}.round1.original_numt.fa"
        ln -sf "\${numt_fai}" "\${SAMPLE_ID}.round1.original_numt.fa.fai"
    fi

    ########################################
    # 3) Recursively find exact round1 split VCF
    #    If multiple Cromwell runs exist, select the newest file by modification time.
    ########################################
    expected_vcf_name="\${SAMPLE_ID}${params.round1_vcf_suffix}"

    mapfile -t vcf_hits < <(
        find "\${VCF_ROOT}" \
            -type f \
            -name "\${expected_vcf_name}" \
            | sort
    )

    if (( \${#vcf_hits[@]} == 0 )); then
        echo "ERROR: Cannot find round1 split VCF for sample \${SAMPLE_ID}" >&2
        echo "Searched recursively under: \${VCF_ROOT}" >&2
        echo "Expected file name: \${expected_vcf_name}" >&2
        exit 1
    fi

    if (( \${#vcf_hits[@]} > 1 )); then
        echo "[WARN] Found multiple round1 split VCF files for sample \${SAMPLE_ID}. Selecting the newest one by modification time:" >&2
        printf '  %s\\n' "\${vcf_hits[@]}" >&2

        vcf="\$(
            find "\${VCF_ROOT}" \
                -type f \
                -name "\${expected_vcf_name}" \
                -printf '%T@\\t%p\\n' \
            | sort -k1,1nr \
            | head -n 1 \
            | cut -f2-
        )"

        echo "[INFO] Selected newest round1 VCF: \${vcf}" >&2
    else
        vcf="\${vcf_hits[0]}"
    fi

    [[ -s "\${vcf}" ]] || {
        echo "ERROR: Selected VCF is missing or empty: \${vcf}" >&2
        exit 1
    }

    echo "[INFO] Selected round1 candidate BAM: \${bam}"
    echo "[INFO] Selected round1 VCF: \${vcf}"

    ########################################
    # 4) Standardize VCF as bgzipped/indexed VCF
    ########################################
    bgzip -c "\${vcf}" > "\${SAMPLE_ID}.round1.input.vcf.gz"
    tabix -p vcf "\${SAMPLE_ID}.round1.input.vcf.gz"

    ########################################
    # 5) Find and standardize Round 1 decoy-coordinate NUMT VCF for consensus NUMT.
    #    This is produced by CALL_NUMT_VARIANTS_DECOY under round_1/numt_decoy_variant_calling.
    #    If older Round 1 outputs do not have it, keep original NUMTs by creating a valid empty VCF.
    ########################################
    mapfile -t nuc_vcf_hits < <(
        if [[ -d "\${NUMT_VCF_ROOT}" ]]; then
            find "\${NUMT_VCF_ROOT}" \
                -type f \
                \\( -name "\${SAMPLE_ID}${params.round1_nuc_vcf_suffix}" -o -name "*${params.round1_nuc_vcf_suffix}" \\) \
                | sort
        fi
    )

    if (( \${#nuc_vcf_hits[@]} == 0 )); then
        echo "[WARN] Cannot find Round 1 decoy NUMT VCF for sample \${SAMPLE_ID} under \${NUMT_VCF_ROOT}; consensus NUMT will equal original NUMT" >&2
        {
          echo '##fileformat=VCFv4.2'
          echo -e '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'
        } | bgzip -c > "\${SAMPLE_ID}.round1.nuc.input.vcf.gz"
        tabix -p vcf -f "\${SAMPLE_ID}.round1.nuc.input.vcf.gz"
    else
        if (( \${#nuc_vcf_hits[@]} > 1 )); then
            echo "[WARN] Found multiple Round 1 decoy NUMT VCF files for sample \${SAMPLE_ID}. Selecting the newest one by modification time:" >&2
            printf '  %s\n' "\${nuc_vcf_hits[@]}" >&2
            nuc_vcf="\$(
                find "\${NUMT_VCF_ROOT}" \
                    -type f \
                    \\( -name "\${SAMPLE_ID}${params.round1_nuc_vcf_suffix}" -o -name "*${params.round1_nuc_vcf_suffix}" \\) \
                    -printf '%T@\t%p\n' \
                | sort -k1,1nr \
                | head -n 1 \
                | cut -f2-
            )"
        else
            nuc_vcf="\${nuc_vcf_hits[0]}"
        fi
        echo "[INFO] Selected round1 decoy NUMT VCF: \${nuc_vcf}"
        if [[ "\${nuc_vcf}" == *.gz ]]; then
            cp "\${nuc_vcf}" "\${SAMPLE_ID}.round1.nuc.input.vcf.gz"
            if [[ -s "\${nuc_vcf}.tbi" ]]; then
                cp "\${nuc_vcf}.tbi" "\${SAMPLE_ID}.round1.nuc.input.vcf.gz.tbi"
            elif [[ -s "\${nuc_vcf%.gz}.idx" ]]; then
                cp "\${nuc_vcf%.gz}.idx" "\${SAMPLE_ID}.round1.nuc.input.vcf.gz.tbi"
            else
                tabix -p vcf -f "\${SAMPLE_ID}.round1.nuc.input.vcf.gz"
            fi
        else
            bgzip -c "\${nuc_vcf}" > "\${SAMPLE_ID}.round1.nuc.input.vcf.gz"
            tabix -p vcf -f "\${SAMPLE_ID}.round1.nuc.input.vcf.gz"
        fi
    fi

    ########################################
    # 6) Sanity checks
    ########################################
    samtools quickcheck -v "\${SAMPLE_ID}.round1.input.bam"
    bcftools index -n "\${SAMPLE_ID}.round1.input.vcf.gz" >/dev/null
    bcftools index -n "\${SAMPLE_ID}.round1.nuc.input.vcf.gz" >/dev/null
    if [[ -s "\${SAMPLE_ID}.round1.original_numt.fa" ]]; then
        samtools faidx "\${SAMPLE_ID}.round1.original_numt.fa" >/dev/null
    else
        echo "[INFO] Skipping faidx sanity check for empty original NUMT FASTA" >&2
    fi

    echo "[INFO] Selected round1 original NUMT FASTA: \${numt_fa}"
    echo "[INFO] Round1 input files prepared for \${SAMPLE_ID}"
    """
}

process BUILD_CONSENSUS_REFERENCE {
    tag { "Build consensus WDL ref bundle for ${meta.id}" }
    label 'alignment_related'
    publishDir "${params.outdir}/${meta.id}/round_2/round2_consensus_ref", mode: 'copy'

    input:
    tuple val(meta), val(species_name), val(ref_name), path(round1_vcf), path(original_numt_fa), path(original_numt_fai), path(round1_nuc_vcf), path(round1_nuc_vcf_tbi)

    output:
    tuple val(meta), val(species_name), val(ref_name),
          path("${meta.id}.round2.consensus.fa"),
          path("${meta.id}.round2.consensus.fa.fai"),
          path("${meta.id}.round2.consensus.dict"),
          path("${meta.id}.round2.consensus.fa.amb"),
          path("${meta.id}.round2.consensus.fa.ann"),
          path("${meta.id}.round2.consensus.fa.bwt"),
          path("${meta.id}.round2.consensus.fa.pac"),
          path("${meta.id}.round2.consensus.fa.sa"),
          path("${meta.id}.round2.consensus.shifted.fa"),
          path("${meta.id}.round2.consensus.shifted.fa.fai"),
          path("${meta.id}.round2.consensus.shifted.dict"),
          path("${meta.id}.round2.consensus.shifted.fa.amb"),
          path("${meta.id}.round2.consensus.shifted.fa.ann"),
          path("${meta.id}.round2.consensus.shifted.fa.bwt"),
          path("${meta.id}.round2.consensus.shifted.fa.pac"),
          path("${meta.id}.round2.consensus.shifted.fa.sa"),
          path("${meta.id}.round2.selfref.fa"),
          path("${meta.id}.round2.selfref.fa.fai"),
          path("${meta.id}.round2.selfref.dict"),
          path("${meta.id}.round2.selfref.fa.amb"),
          path("${meta.id}.round2.selfref.fa.ann"),
          path("${meta.id}.round2.selfref.fa.bwt"),
          path("${meta.id}.round2.selfref.fa.pac"),
          path("${meta.id}.round2.selfref.fa.sa"),
          path("${meta.id}.round2.selfref.shifted.fa"),
          path("${meta.id}.round2.selfref.shifted.fa.fai"),
          path("${meta.id}.round2.selfref.shifted.dict"),
          path("${meta.id}.round2.selfref.shifted.fa.amb"),
          path("${meta.id}.round2.selfref.shifted.fa.ann"),
          path("${meta.id}.round2.selfref.shifted.fa.bwt"),
          path("${meta.id}.round2.selfref.shifted.fa.pac"),
          path("${meta.id}.round2.selfref.shifted.fa.sa"),
          path("${meta.id}.round2.consensus.non_control_region.interval_list"),
          path("${meta.id}.round2.consensus.control_region_shifted.interval_list"),
          path("${meta.id}.round2.consensus.ShiftBack.chain"),
          path("${meta.id}.round2.consensus_sites.vcf.gz"),
          path("${meta.id}.round2.consensus_sites.vcf.gz.tbi"),
          emit: consensus_ref

    script:
    def mt_fasta = "${params.ref_dir}/${ref_name}.fasta"
    def numt_fasta = original_numt_fa
    def strictNumtRef = params.strict_numt_ref ? "true" : "false"
    def filterExpr = params.consensus_filter_expr.toString().replace("'", "'\\''")
    """
    #!/usr/bin/env bash
    set -euo pipefail

    MT_FASTA="${mt_fasta}"
    FILTER_EXPR='${filterExpr}'
    PICARD="${params.picard_jar}"
    SHIFT=${params.mt_shift}
    PADDING=${params.mt_interval_padding}
    LEFT_MARGIN=${params.mt_left_margin}
    RIGHT_MARGIN=${params.mt_right_margin}

    [[ -s "\${MT_FASTA}" ]] || { echo "ERROR: Missing mt FASTA: \${MT_FASTA}" >&2; exit 1; }
    [[ -s "\${PICARD}" ]] || { echo "ERROR: Missing Picard jar: \${PICARD}" >&2; exit 1; }

    index_fasta() {
        local fa="\$1"
        local dict="\$2"
        local fai="\${fa}.fai"
        local bwa_files=("\${fa}.amb" "\${fa}.ann" "\${fa}.bwt" "\${fa}.pac" "\${fa}.sa")

        rm -f "\${fai}" "\${dict}" "\${bwa_files[@]}"
        echo "[INFO] samtools faidx: \${fa}"
        samtools faidx "\${fa}"
        echo "[INFO] bwa index: \${fa}"
        bwa index "\${fa}"
        echo "[INFO] Picard CreateSequenceDictionary: \${fa}"
        java -jar "\${PICARD}" CreateSequenceDictionary \
            R="\${fa}" \
            O="\${dict}"
    }

    # 1) Normalize and filter round1 VCF.
    #    Use high-confidence near-homoplasmic SNVs and indels for consensus construction.
    if [[ "${round1_vcf}" == *.vcf.gz ]]; then
        ln -sf ${round1_vcf} ${meta.id}.round1.vcf.gz
        bcftools index -t -f ${meta.id}.round1.vcf.gz
    else
        bgzip -c ${round1_vcf} > ${meta.id}.round1.vcf.gz
        tabix -p vcf ${meta.id}.round1.vcf.gz
    fi

    bcftools norm -m-any -f "\${MT_FASTA}" ${meta.id}.round1.vcf.gz -Oz -o ${meta.id}.round1.norm.vcf.gz
    bcftools index -t -f ${meta.id}.round1.norm.vcf.gz

    bcftools view \
        -f PASS \
        -m2 -M2 \
        -v snps,indels \
        -i "\${FILTER_EXPR}" \
        ${meta.id}.round1.norm.vcf.gz \
        -Oz -o ${meta.id}.round2.consensus_sites.vcf.gz
    bcftools index -t -f ${meta.id}.round2.consensus_sites.vcf.gz
    bcftools view -H ${meta.id}.round2.consensus_sites.vcf.gz | wc -l > ${meta.id}.round2.consensus_sites.n.txt

    # 2) Build unshifted consensus FASTA from round1 SNVs and indels.
    bcftools consensus -H A -f "\${MT_FASTA}" ${meta.id}.round2.consensus_sites.vcf.gz > ${meta.id}.round2.consensus.fa

    # 3) Generate shifted FASTA, interval lists, and ShiftBack.chain directly inside NF.
    #    Coordinate convention follows the chain/interval logic:
    #      shifted FASTA = original bases (SHIFT+1..end) + (1..SHIFT)
    #      ShiftBack.chain maps shifted-reference coordinates back to unshifted coordinates.
    python3 - <<'PY_IN_NF'
from pathlib import Path
import hashlib

sample = "${meta.id}"
fa = Path(f"{sample}.round2.consensus.fa")
shifted_fa = Path(f"{sample}.round2.consensus.shifted.fa")
non_control = Path(f"{sample}.round2.consensus.non_control_region.interval_list")
control_shifted = Path(f"{sample}.round2.consensus.control_region_shifted.interval_list")
chain = Path(f"{sample}.round2.consensus.ShiftBack.chain")
shift = int("${params.mt_shift}")
padding = int("${params.mt_interval_padding}")
left_margin = int("${params.mt_left_margin}")
right_margin = int("${params.mt_right_margin}")

name = None
seq_parts = []
for line in fa.read_text().splitlines():
    if line.startswith(">"):
        if name is not None:
            raise SystemExit("ERROR: consensus FASTA must contain exactly one contig")
        name = line[1:].split()[0]
    else:
        seq_parts.append(line.strip())

if name is None:
    raise SystemExit("ERROR: empty consensus FASTA")
seq = "".join(seq_parts).upper()
L = len(seq)
if L <= 0:
    raise SystemExit("ERROR: consensus sequence length is zero")
if not (0 < shift < L):
    raise SystemExit(f"ERROR: mt_shift {shift} must be between 1 and length-1 (length={L})")

# B/official WDL shift: original coordinate SHIFT+1 becomes shifted coordinate 1.
shifted = seq[shift:] + seq[:shift]
with shifted_fa.open("w") as handle:
    handle.write(f">{name}\\n")
    for i in range(0, len(shifted), 60):
        handle.write(shifted[i:i+60] + "\\n")

def md5_seq(s: str) -> str:
    return hashlib.md5(s.upper().encode()).hexdigest()

def abs_uri(p: Path) -> str:
    return "file://" + str(p.resolve())

nc_start = left_margin + 1
nc_end = L - right_margin
if nc_start < 1 or nc_end > L or nc_start > nc_end:
    raise SystemExit(f"ERROR: invalid non-control interval: start={nc_start}, end={nc_end}, length={L}")

if padding >= nc_start or (nc_end + padding) > L:
    print(f"[WARN] Padding {padding} is too large for length {L}; using padding=0")
    P = 0
else:
    P = padding

# Same mapping as make_mt_intervals.sh shift_forward():
# original coordinate x -> shifted coordinate ((x - 1 - SHIFT) mod L) + 1
# for the left-shifted reference generated above.
def shift_forward(x: int) -> int:
    return ((x - 1 - shift) % L) + 1

with non_control.open("w") as handle:
    handle.write("@HD\\tVN:1.6\\n")
    handle.write(f"@SQ\\tSN:{name}\\tLN:{L}\\tM5:{md5_seq(seq)}\\tUR:{abs_uri(fa)}\\n")
    handle.write(f"{name}\\t{nc_start - P}\\t{nc_end + P}\\t+\\t.\\n")

start_raw = nc_end + 1 - P
end_raw = nc_start - 1 + P
start_shifted = shift_forward(start_raw)
end_shifted = shift_forward(end_raw)
with control_shifted.open("w") as handle:
    handle.write("@HD\\tVN:1.6\\n")
    handle.write(f"@SQ\\tSN:{name}\\tLN:{L}\\tM5:{md5_seq(shifted)}\\tUR:{abs_uri(shifted_fa)}\\n")
    handle.write(f"{name}\\t{start_shifted}\\t{end_shifted}\\t+\\t.\\n")

left = L - shift
with chain.open("w") as handle:
    handle.write(f"chain {left} {name} {L} + 0 {left} {name} {L} + {shift} {L} 1\\n")
    handle.write(f"{left}\\n\\n")
    handle.write(f"chain {shift} {name} {L} + {left} {L} {name} {L} + 0 {shift} 2\\n")
    handle.write(f"{shift}\\n")

print(f"[INFO] Wrote shifted FASTA: {shifted_fa}")
print(f"[INFO] Wrote interval lists: {non_control}, {control_shifted} (padding={P})")
print(f"[INFO] Wrote ShiftBack chain: {chain}")
PY_IN_NF

    # 4) Build consensus NUMT FASTA and mtSwirl-like self-reference FASTAs.
    #    NUMT contigs start from Round 1 sample-specific original_numt.fa and are updated
    #    with the Round 1 NUMT HaplotypeCaller VCF. If no NUMT VCF records are available,
    #    the consensus NUMT FASTA is identical to the original NUMT FASTA.
    python3 - <<'PY_SELFREF'
from pathlib import Path
import gzip
import sys

sample = "${meta.id}"
numt_path = Path("${numt_fasta}") if "${numt_fasta}" else None
nuc_vcf_path = Path("${round1_nuc_vcf}")
strict = "${strictNumtRef}".lower() == "true"

standard_chrM = Path(f"{sample}.round2.consensus.fa")
shifted_chrM = Path(f"{sample}.round2.consensus.shifted.fa")
consensus_numt = Path(f"{sample}.round2.consensus_numt.fa")
standard_selfref = Path(f"{sample}.round2.selfref.fa")
shifted_selfref = Path(f"{sample}.round2.selfref.shifted.fa")

# Parse first contig name from consensus chrM so NUMT headers can be made unique.
chr_name = None
for line in standard_chrM.read_text().splitlines():
    if line.startswith(">"): 
        chr_name = line[1:].split()[0]
        break
if chr_name is None:
    raise SystemExit("ERROR: cannot determine chrM contig name from consensus FASTA")

seen = {chr_name}

def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def write_wrapped(handle, seq, width=60):
    seq = seq.replace(" ", "").replace("\\t", "").upper()
    for i in range(0, len(seq), width):
        handle.write(seq[i:i+width] + "\\n")


def read_fasta_records(path):
    records = []
    header = None
    seq_parts = []
    with Path(path).open() as handle:
        for line in handle:
            line = line.rstrip("\\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            records.append((header, "".join(seq_parts).upper()))
    return records


def parse_numt_header(header):
    # Round 1 PREPARE_DECOY_REFERENCE names intervals as NUMT_<n>_<chrom>_<1-based-start>_<end>.
    raw_name = header.split()[0]
    parts = raw_name.split("_")
    if len(parts) < 5 or parts[0] != "NUMT":
        return None
    try:
        start = int(parts[-2])
        end = int(parts[-1])
    except ValueError:
        return None
    chrom = "_".join(parts[2:-2])
    if not chrom or start < 1 or end < start:
        return None
    return chrom, start, end


def read_vcf_records(path):
    by_chrom = {}
    if not Path(path).exists():
        return by_chrom
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) < 8:
                continue
            chrom, pos_s, _id, ref, alt, qual, filt, info = fields[:8]
            if filt not in ("PASS", "."):
                continue
            if "," in alt or alt in (".", "*") or alt.startswith("<"):
                continue
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            by_chrom.setdefault(chrom, []).append((pos, ref.upper(), alt.upper()))
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda x: x[0])
    return by_chrom


def build_consensus_numt(original_fa, nuc_vcf, out_fa):
    if original_fa is None or str(original_fa) == "":
        print("[WARN] original NUMT FASTA was not provided; consensus NUMT FASTA is empty", file=sys.stderr)
        out_fa.write_text("")
        return 0, 0
    if not original_fa.exists():
        msg = f"NUMT FASTA not found: {original_fa}"
        if strict:
            raise SystemExit("ERROR: " + msg)
        print("[WARN] " + msg + "; consensus NUMT FASTA is empty", file=sys.stderr)
        out_fa.write_text("")
        return 0, 0
    if original_fa.stat().st_size == 0:
        print(f"[INFO] NUMT FASTA is empty: {original_fa}; self-reference is chrM-only for this sample", file=sys.stderr)
        out_fa.write_text("")
        return 0, 0

    variants = read_vcf_records(nuc_vcf)
    records = read_fasta_records(original_fa)
    applied_total = 0

    with out_fa.open("w") as out:
        for header, seq in records:
            parsed = parse_numt_header(header)
            seq_list = list(seq)
            applied = 0
            offset = 0
            last_ref_end = -1
            if parsed is not None:
                chrom, start, end = parsed
                for pos, ref, alt in variants.get(chrom, []):
                    if pos < start or pos > end:
                        continue
                    ref_end = pos + len(ref) - 1
                    if ref_end > end:
                        continue
                    if pos <= last_ref_end:
                        print(f"[WARN] Skipping overlapping NUMT variant {chrom}:{pos} {ref}>{alt} for {header}", file=sys.stderr)
                        continue
                    local0 = pos - start + offset
                    observed = "".join(seq_list[local0:local0 + len(ref)]).upper()
                    if observed != ref:
                        print(
                            f"[WARN] Skipping NUMT variant {chrom}:{pos} {ref}>{alt} for {header}: observed {observed}",
                            file=sys.stderr,
                        )
                        continue
                    seq_list[local0:local0 + len(ref)] = list(alt)
                    offset += len(alt) - len(ref)
                    last_ref_end = ref_end
                    applied += 1
            else:
                print(f"[WARN] Cannot parse NUMT FASTA header; leaving sequence unchanged: {header}", file=sys.stderr)

            out.write(f">{header}\\n")
            write_wrapped(out, "".join(seq_list))
            applied_total += applied
    print(f"[INFO] Wrote consensus NUMT FASTA: {out_fa} (records={len(records)}, applied_variants={applied_total})")
    return len(records), applied_total


def append_sanitized_numts(out_handle, numt_fa):
    if numt_fa is None or str(numt_fa) == "" or not numt_fa.exists() or numt_fa.stat().st_size == 0:
        print("[INFO] consensus NUMT FASTA is empty; self-reference is chrM-only for this sample", file=sys.stderr)
        return 0

    n = 0
    for header, seq in read_fasta_records(numt_fa):
        raw_name = header.split()[0]
        clean = "".join(c if c.isalnum() or c in "._:-" else "_" for c in raw_name)
        if not clean:
            clean = f"NUMT_{n+1}"
        if clean == chr_name or clean in seen:
            clean = f"NUMT_{n+1}_{clean}"
        while clean in seen:
            clean = f"NUMT_{n+1}_{clean}"
        seen.add(clean)
        if seq:
            out_handle.write(f">{clean} consensus_NUMT source={raw_name}\\n")
            write_wrapped(out_handle, seq)
            n += 1
    print(f"[INFO] Appended {n} consensus NUMT contigs from {numt_fa}")
    return n

build_consensus_numt(numt_path, nuc_vcf_path, consensus_numt)

for chr_fa, self_fa in [(standard_chrM, standard_selfref), (shifted_chrM, shifted_selfref)]:
    seen.clear(); seen.add(chr_name)
    chr_text = chr_fa.read_text()
    with self_fa.open("w") as out:
        out.write(chr_text)
        if not chr_text.endswith("\\n"):
            out.write("\\n")
        append_sanitized_numts(out, consensus_numt)

PY_SELFREF

    # 5) Index chrM-only consensus references and combined self-references.
    index_fasta ${meta.id}.round2.consensus.fa ${meta.id}.round2.consensus.dict
    index_fasta ${meta.id}.round2.consensus.shifted.fa ${meta.id}.round2.consensus.shifted.dict
    index_fasta ${meta.id}.round2.selfref.fa ${meta.id}.round2.selfref.dict
    index_fasta ${meta.id}.round2.selfref.shifted.fa ${meta.id}.round2.selfref.shifted.dict

    # 6) Rebuild interval lists so their sequence dictionaries match the combined self-references.
    #    Intervals still target only chrM, but @SQ lines include chrM + NUMT contigs.
    python3 - <<'PY_INTERVALS'
from pathlib import Path

sample = "${meta.id}"
shift = int("${params.mt_shift}")
padding = int("${params.mt_interval_padding}")
left_margin = int("${params.mt_left_margin}")
right_margin = int("${params.mt_right_margin}")

cons_fai = Path(f"{sample}.round2.consensus.fa.fai")
name, L = None, None
with cons_fai.open() as handle:
    first = handle.readline().rstrip("\\n").split("\\t")
    name, L = first[0], int(first[1])

nc_start = left_margin + 1
nc_end = L - right_margin
if nc_start < 1 or nc_end > L or nc_start > nc_end:
    raise SystemExit(f"ERROR: invalid non-control interval: start={nc_start}, end={nc_end}, length={L}")
P = 0 if (padding >= nc_start or (nc_end + padding) > L) else padding

def shift_forward(x):
    return ((x - 1 - shift) % L) + 1

non_start = nc_start - P
non_end = nc_end + P
start_raw = nc_end + 1 - P
end_raw = nc_start - 1 + P
control_start = shift_forward(start_raw)
control_end = shift_forward(end_raw)

def sq_lines_from_dict(dict_path):
    lines = ["@HD\\tVN:1.6\\n"]
    with Path(dict_path).open() as handle:
        for line in handle:
            if line.startswith("@SQ"):
                lines.append(line)
    return lines

def write_interval(out_path, dict_path, contig, start, end):
    with Path(out_path).open("w") as out:
        out.writelines(sq_lines_from_dict(dict_path))
        out.write(f"{contig}\\t{start}\\t{end}\\t+\\t.\\n")

write_interval(
    f"{sample}.round2.consensus.non_control_region.interval_list",
    f"{sample}.round2.selfref.dict",
    name,
    non_start,
    non_end,
)
write_interval(
    f"{sample}.round2.consensus.control_region_shifted.interval_list",
    f"{sample}.round2.selfref.shifted.dict",
    name,
    control_start,
    control_end,
)
print(f"[INFO] Rewrote interval lists with combined self-reference dictionaries; padding={P}")
PY_INTERVALS

    # 7) Basic sanity checks for all WDL reference files generated here.
    for f in \
        ${meta.id}.round2.consensus.fa \
        ${meta.id}.round2.consensus.fa.fai \
        ${meta.id}.round2.consensus.dict \
        ${meta.id}.round2.consensus.fa.amb \
        ${meta.id}.round2.consensus.fa.ann \
        ${meta.id}.round2.consensus.fa.bwt \
        ${meta.id}.round2.consensus.fa.pac \
        ${meta.id}.round2.consensus.fa.sa \
        ${meta.id}.round2.consensus.shifted.fa \
        ${meta.id}.round2.consensus.shifted.fa.fai \
        ${meta.id}.round2.consensus.shifted.dict \
        ${meta.id}.round2.consensus.shifted.fa.amb \
        ${meta.id}.round2.consensus.shifted.fa.ann \
        ${meta.id}.round2.consensus.shifted.fa.bwt \
        ${meta.id}.round2.consensus.shifted.fa.pac \
        ${meta.id}.round2.consensus.shifted.fa.sa \
        ${meta.id}.round2.selfref.fa \
        ${meta.id}.round2.selfref.fa.fai \
        ${meta.id}.round2.selfref.dict \
        ${meta.id}.round2.selfref.fa.amb \
        ${meta.id}.round2.selfref.fa.ann \
        ${meta.id}.round2.selfref.fa.bwt \
        ${meta.id}.round2.selfref.fa.pac \
        ${meta.id}.round2.selfref.fa.sa \
        ${meta.id}.round2.selfref.shifted.fa \
        ${meta.id}.round2.selfref.shifted.fa.fai \
        ${meta.id}.round2.selfref.shifted.dict \
        ${meta.id}.round2.selfref.shifted.fa.amb \
        ${meta.id}.round2.selfref.shifted.fa.ann \
        ${meta.id}.round2.selfref.shifted.fa.bwt \
        ${meta.id}.round2.selfref.shifted.fa.pac \
        ${meta.id}.round2.selfref.shifted.fa.sa \
        ${meta.id}.round2.consensus.non_control_region.interval_list \
        ${meta.id}.round2.consensus.control_region_shifted.interval_list \
        ${meta.id}.round2.consensus.ShiftBack.chain \
        ${meta.id}.round2.consensus_sites.vcf.gz \
        ${meta.id}.round2.consensus_sites.vcf.gz.tbi
    do
        [[ -s "\$f" ]] || { echo "ERROR: expected WDL reference file missing or empty: \$f" >&2; exit 1; }
    done
    """
}

process REALIGN_TO_CONSENSUS_ASSIGNED_BAMS {
    tag { "Realign reads to standard/shifted consensus and extract chrM-assigned reads for ${meta.id}" }
    label 'alignment_related'
    publishDir "${params.outdir}/${meta.id}/round_2/round2_alignment_consensus", mode: 'copy', pattern: "*.{bam,bai,txt,log}"

    input:
    tuple val(meta), val(species_name), val(ref_name), path(round1_bam),
          path(consensus_fa), path(consensus_fai), path(consensus_dict),
          path(consensus_amb), path(consensus_ann), path(consensus_bwt), path(consensus_pac), path(consensus_sa),
          path(shifted_fa), path(shifted_fai), path(shifted_dict),
          path(shifted_amb), path(shifted_ann), path(shifted_bwt), path(shifted_pac), path(shifted_sa),
          path(selfref_fa), path(selfref_fai), path(selfref_dict),
          path(selfref_amb), path(selfref_ann), path(selfref_bwt), path(selfref_pac), path(selfref_sa),
          path(selfref_shifted_fa), path(selfref_shifted_fai), path(selfref_shifted_dict),
          path(selfref_shifted_amb), path(selfref_shifted_ann), path(selfref_shifted_bwt), path(selfref_shifted_pac), path(selfref_shifted_sa),
          path(non_control_interval), path(control_shifted_interval), path(shift_back_chain),
          path(consensus_sites_vcf), path(consensus_sites_tbi)

    output:
    tuple val(meta), val(species_name), val(ref_name),
          path("${meta.id}.standard.chrM_assigned.bam"),
          path("${meta.id}.standard.chrM_assigned.bam.bai"),
          path("${meta.id}.shifted.chrM_assigned.bam"),
          path("${meta.id}.shifted.chrM_assigned.bam.bai"),
          emit: assigned_bams

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    THREADS=${task.cpus ?: 4}
    SAMPLE_ID="${meta.id}"
    PICARD="${params.picard_jar}"
    CHR_NAME=\$(awk 'NR==1{print \$1}' ${consensus_fai})

    samtools quickcheck -v ${round1_bam}

    echo "[INFO] Extract FASTQ from round1 candidate BAM: ${round1_bam}"
    samtools collate -@ "\${THREADS}" -u -O ${round1_bam} "\${SAMPLE_ID}.collate.tmp" \
      | samtools fastq -@ "\${THREADS}" -n -F 0x900 \
          -1 "\${SAMPLE_ID}.R1.fq" \
          -2 "\${SAMPLE_ID}.R2.fq" \
          -0 /dev/null \
          -s "\${SAMPLE_ID}.singleton.fq" \
          -

    RG="@RG\\tID:\${SAMPLE_ID}\\tSM:\${SAMPLE_ID}\\tPL:ILLUMINA\\tLB:\${SAMPLE_ID}"

    align_branch() {
        local branch="\$1"
        local ref_fa="\$2"
        local chr_name="\$3"

        echo "[INFO] Align \${branch} branch to \${ref_fa}"

        bwa mem -K 100000000 -v 3 -t "\${THREADS}" -Y -R "\${RG}" \
            "\${ref_fa}" \
            "\${SAMPLE_ID}.R1.fq" "\${SAMPLE_ID}.R2.fq" \
            2> "\${SAMPLE_ID}.\${branch}.bwa.paired.stderr.log" \
          | samtools sort -@ "\${THREADS}" \
              -o "\${SAMPLE_ID}.\${branch}.selfref.paired.sorted.bam" -

        if [[ -s "\${SAMPLE_ID}.singleton.fq" ]] && [[ \$(wc -l < "\${SAMPLE_ID}.singleton.fq") -gt 0 ]]; then
            bwa mem -K 100000000 -v 3 -t "\${THREADS}" -Y -R "\${RG}" \
                "\${ref_fa}" \
                "\${SAMPLE_ID}.singleton.fq" \
                2> "\${SAMPLE_ID}.\${branch}.bwa.singleton.stderr.log" \
              | samtools sort -@ "\${THREADS}" \
                  -o "\${SAMPLE_ID}.\${branch}.selfref.singleton.sorted.bam" -

            samtools merge -@ "\${THREADS}" -f \
                "\${SAMPLE_ID}.\${branch}.selfref.sorted.bam" \
                "\${SAMPLE_ID}.\${branch}.selfref.paired.sorted.bam" \
                "\${SAMPLE_ID}.\${branch}.selfref.singleton.sorted.bam"
        else
            mv "\${SAMPLE_ID}.\${branch}.selfref.paired.sorted.bam" \
               "\${SAMPLE_ID}.\${branch}.selfref.sorted.bam"
        fi

        java -Xmx8G -jar "\${PICARD}" MarkDuplicates \
            INPUT="\${SAMPLE_ID}.\${branch}.selfref.sorted.bam" \
            OUTPUT="\${SAMPLE_ID}.\${branch}.selfref.md.bam" \
            METRICS_FILE="\${SAMPLE_ID}.\${branch}.duplicate_metrics.txt" \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=LENIENT \
            ASSUME_SORT_ORDER=coordinate

        # mtSwirl-like: after competitive mapping/preprocessing, keep reads mapping to chrM.
        # Do NOT require proper pair or mate-on-same-contig here.
        samtools view -@ "\${THREADS}" -b -F 2308 \
            "\${SAMPLE_ID}.\${branch}.selfref.md.bam" \
            "\${chr_name}" \
          > "\${SAMPLE_ID}.\${branch}.chrM_assigned.bam"

        samtools index -@ "\${THREADS}" "\${SAMPLE_ID}.\${branch}.chrM_assigned.bam"
        samtools quickcheck -v "\${SAMPLE_ID}.\${branch}.chrM_assigned.bam"
    }

    align_branch "standard" "${selfref_fa}" "\${CHR_NAME}"
    align_branch "shifted"  "${selfref_shifted_fa}" "\${CHR_NAME}"

    echo "[INFO] Finished mtSwirl-like preassigned BAMs for \${SAMPLE_ID}"
    """
}

process GENERATE_BAM_TSV {
    tag { "Generate 4-column assigned-BAM TSV for ${meta.id}" }
    label 'generation_related'
    publishDir "${params.outdir}/${meta.id}/round_2/round2_wdl_inputs", mode: 'copy'

    input:
    tuple val(meta), val(species_name), val(ref_name),
          path(std_bam), path(std_bai), path(shift_bam), path(shift_bai)

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}_round2_assigned_bam_list.tsv"), emit: tsv

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    echo -e "\$(readlink -f ${std_bam})\t\$(readlink -f ${std_bai})\t\$(readlink -f ${shift_bam})\t\$(readlink -f ${shift_bai})" > ${meta.id}_round2_assigned_bam_list.tsv
    """
}

process GENERATE_WDL_JSON_ROUND2 {
    tag { "Generate round2 WDL JSON for ${meta.id}" }
    label 'generation_related'
    publishDir "${params.outdir}/${meta.id}/round_2/round2_wdl_inputs", mode: 'copy'

    input:
    tuple val(meta), val(species_name), val(ref_name), path(cram_tsv),
          path(consensus_fa), path(consensus_fai), path(consensus_dict),
          path(consensus_amb), path(consensus_ann), path(consensus_bwt), path(consensus_pac), path(consensus_sa),
          path(shifted_fa), path(shifted_fai), path(shifted_dict),
          path(shifted_amb), path(shifted_ann), path(shifted_bwt), path(shifted_pac), path(shifted_sa),
          path(selfref_fa), path(selfref_fai), path(selfref_dict),
          path(selfref_amb), path(selfref_ann), path(selfref_bwt), path(selfref_pac), path(selfref_sa),
          path(selfref_shifted_fa), path(selfref_shifted_fai), path(selfref_shifted_dict),
          path(selfref_shifted_amb), path(selfref_shifted_ann), path(selfref_shifted_bwt), path(selfref_shifted_pac), path(selfref_shifted_sa),
          path(non_control_interval), path(control_shifted_interval), path(shift_back_chain),
          path(consensus_sites_vcf), path(consensus_sites_tbi)

    output:
    tuple val(meta), val(species_name), val(ref_name), path("${meta.id}_round2_wdl_inputs.json"), emit: json

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    CHR_NAME=\$(awk 'NR==1{print \$1}' ${consensus_fai})
    MT_LEN=\$(awk 'NR==1{print \$2}' ${consensus_fai})

    cat > ${meta.id}_round2_wdl_inputs.json <<EOFJSON
{
  "MitochondriaMultiSamplePipeline.picard": "${params.picard_jar}",
  "MitochondriaMultiSamplePipeline.haplocheckCLI": "${params.haplocheck_jar}",
  "MitochondriaMultiSamplePipeline.gatk": "${params.gatk_jar}",
  "MitochondriaMultiSamplePipeline.compress_output_vcf": true,
  "MitochondriaMultiSamplePipeline.inputSamplesFile": "\$(readlink -f ${cram_tsv})",

  "MitochondriaMultiSamplePipeline.ref_fasta": "\$(readlink -f ${selfref_fa})",
  "MitochondriaMultiSamplePipeline.ref_fasta_index": "\$(readlink -f ${selfref_fai})",
  "MitochondriaMultiSamplePipeline.ref_dict": "\$(readlink -f ${selfref_dict})",

  "MitochondriaMultiSamplePipeline.mt_dict": "\$(readlink -f ${selfref_dict})",
  "MitochondriaMultiSamplePipeline.mt_fasta": "\$(readlink -f ${selfref_fa})",
  "MitochondriaMultiSamplePipeline.mt_fasta_index": "\$(readlink -f ${selfref_fai})",
  "MitochondriaMultiSamplePipeline.mt_amb": "\$(readlink -f ${selfref_amb})",
  "MitochondriaMultiSamplePipeline.mt_ann": "\$(readlink -f ${selfref_ann})",
  "MitochondriaMultiSamplePipeline.mt_bwt": "\$(readlink -f ${selfref_bwt})",
  "MitochondriaMultiSamplePipeline.mt_pac": "\$(readlink -f ${selfref_pac})",
  "MitochondriaMultiSamplePipeline.mt_sa": "\$(readlink -f ${selfref_sa})",

  "MitochondriaMultiSamplePipeline.mt_shifted_dict": "\$(readlink -f ${selfref_shifted_dict})",
  "MitochondriaMultiSamplePipeline.mt_shifted_fasta": "\$(readlink -f ${selfref_shifted_fa})",
  "MitochondriaMultiSamplePipeline.mt_shifted_fasta_index": "\$(readlink -f ${selfref_shifted_fai})",
  "MitochondriaMultiSamplePipeline.mt_shifted_amb": "\$(readlink -f ${selfref_shifted_amb})",
  "MitochondriaMultiSamplePipeline.mt_shifted_ann": "\$(readlink -f ${selfref_shifted_ann})",
  "MitochondriaMultiSamplePipeline.mt_shifted_bwt": "\$(readlink -f ${selfref_shifted_bwt})",
  "MitochondriaMultiSamplePipeline.mt_shifted_pac": "\$(readlink -f ${selfref_shifted_pac})",
  "MitochondriaMultiSamplePipeline.mt_shifted_sa": "\$(readlink -f ${selfref_shifted_sa})",

  "MitochondriaMultiSamplePipeline.shift_back_chain": "\$(readlink -f ${shift_back_chain})",
  "MitochondriaMultiSamplePipeline.non_control_region_interval_list": "\$(readlink -f ${non_control_interval})",
  "MitochondriaMultiSamplePipeline.control_region_shifted_reference_interval_list": "\$(readlink -f ${control_shifted_interval})",

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
    tag "Round2 Variant Calling on ${meta.id}"
    label 'wdl_related'
    publishDir "${params.outdir}/${meta.id}/round_2/round_2_variant_calling_consensus/", mode: 'copy'

    input:
    tuple val(meta), val(species_name), val(ref_name), path(wdl_inputs_json)
    path cromwell_config

    output:
    tuple val(meta), val(species_name), val(ref_name), path("cromwell-executions"), emit: cromwell_out

    script:
    """
    #!/bin/bash
    set -euo pipefail

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

    java -Dconfig.file=${cromwell_config} \
      -jar ${params.cromwell_jar} run \
      ${params.wdl_script_round2} \
      --inputs ${wdl_inputs_json} \
      --options cromwell_options.json
    """
}

process LIFTBACK_ROUND2_VCF_TO_ORIGINAL {
    tag { "Lift back round2 VCF to original mt coordinates for ${meta.id}" }
    label 'generation_related'
    publishDir "${params.outdir}/${meta.id}/round_2_variant_calling_original_coords", mode: 'copy'

    input:
    tuple val(meta), val(species_name), val(ref_name), path(cromwell_dir),
          path(consensus_fa), path(consensus_fai), path(consensus_dict),
          path(consensus_amb), path(consensus_ann), path(consensus_bwt), path(consensus_pac), path(consensus_sa),
          path(shifted_fa), path(shifted_fai), path(shifted_dict),
          path(shifted_amb), path(shifted_ann), path(shifted_bwt), path(shifted_pac), path(shifted_sa),
          path(selfref_fa), path(selfref_fai), path(selfref_dict),
          path(selfref_amb), path(selfref_ann), path(selfref_bwt), path(selfref_pac), path(selfref_sa),
          path(selfref_shifted_fa), path(selfref_shifted_fai), path(selfref_shifted_dict),
          path(selfref_shifted_amb), path(selfref_shifted_ann), path(selfref_shifted_bwt), path(selfref_shifted_pac), path(selfref_shifted_sa),
          path(non_control_interval), path(control_shifted_interval), path(shift_back_chain),
          path(consensus_sites_vcf), path(consensus_sites_tbi)

    output:
    tuple val(meta), val(species_name), val(ref_name),
          path("${meta.id}.round2.original_coords.clean.final.split.vcf.gz"),
          path("${meta.id}.round2.original_coords.clean.final.split.vcf.gz.tbi"),
          path("${meta.id}.round2.original_coords.liftback.summary.tsv"),
          path("${meta.id}.round2.original_coords.unresolved.tsv"),
          path("${meta.id}.round2.original_coords.skipped_overlap_consensus_sites.tsv"),
          path("${meta.id}.round2.original_coords.round2_dropped_by_round1.tsv"),
          path("${meta.id}.round2.original_coords.merge_policy.summary.tsv"),
          emit: original_coord_vcf

    tuple val(meta), val(species_name), val(ref_name),
          path("${meta.id}.round2.consensus_coords.per_base_coverage.tsv"),
          path("${meta.id}.round2.original_coords.per_base_coverage.tsv"),
          path("${meta.id}.round2.original_coords.inserted_base_coverage.tsv"),
          path("${meta.id}.round2.original_coords.deleted_positions_without_consensus_coverage.tsv"),
          path("${meta.id}.round2.original_coords.coverage_liftback.summary.tsv"),
          emit: lifted_coverage

    script:
    def mt_fasta = "${params.ref_dir}/${ref_name}.fasta"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    SAMPLE_ID="${meta.id}"
    MT_FASTA="${mt_fasta}"
    FINAL_SUFFIX="${params.round2_final_vcf_suffix}"

    [[ -d ${cromwell_dir} ]] || { echo "ERROR: Missing Cromwell directory: ${cromwell_dir}" >&2; exit 1; }
    [[ -s "\${MT_FASTA}" ]] || { echo "ERROR: Missing original mt FASTA: \${MT_FASTA}" >&2; exit 1; }
    [[ -s ${consensus_sites_vcf} ]] || { echo "ERROR: Missing consensus-sites VCF: ${consensus_sites_vcf}" >&2; exit 1; }

    # Find the final split VCF produced by the WDL. Prefer the exact sample suffix,
    # but fall back to any clean.final.split.vcf(.gz) if the WDL output path changes.
    # Use Python instead of a complex find expression because Nextflow/Groovy can
    # mis-parse escaped shell parentheses during script compilation.
    python3 - <<'PY_FIND_VCF' > round2_vcf_path.txt
from pathlib import Path
import sys

sample = "${meta.id}"
suffix = "${params.round2_final_vcf_suffix}"
root = Path("${cromwell_dir}")

if not root.exists():
    print(f"ERROR: Missing Cromwell directory: {root}", file=sys.stderr)
    sys.exit(1)

candidates = []
for p in root.rglob("*"):
    if not p.is_file():
        continue
    name = p.name
    if name.endswith((".idx", ".tbi")):
        continue
    if (
        name == sample + suffix
        or name == sample + suffix + ".gz"
        or name.endswith("clean.final.split.vcf")
        or name.endswith("clean.final.split.vcf.gz")
        or name.endswith("final.split.vcf")
        or name.endswith("final.split.vcf.gz")
    ):
        candidates.append(p)

if not candidates:
    print(f"ERROR: Cannot find final round2 split VCF under {root}", file=sys.stderr)
    print(f"Expected suffix: {suffix}", file=sys.stderr)
    for p in sorted(x for x in root.rglob("*") if x.is_file()):
        print(f"[CROMWELL_FILE] {p}", file=sys.stderr)
    sys.exit(1)

candidates.sort(key=lambda x: x.stat().st_mtime, reverse=True)
if len(candidates) > 1:
    print("[WARN] Found multiple candidate final VCFs. Selecting newest:", file=sys.stderr)
    for p in candidates:
        print(f"  {p}", file=sys.stderr)

print(candidates[0])
PY_FIND_VCF

    ROUND2_VCF="\$(cat round2_vcf_path.txt)"
    echo "[INFO] Selected round2 consensus-coordinate VCF: \${ROUND2_VCF}"
    echo "[INFO] Original mt FASTA: \${MT_FASTA}"
    echo "[INFO] Consensus-sites VCF used for lift-back: ${consensus_sites_vcf}"

    # Work with a local copy of the original mt FASTA so faidx can be created
    # even when the reference directory is read-only.
    cp "\${MT_FASTA}" "\${SAMPLE_ID}.original_mt.fa"
    samtools faidx "\${SAMPLE_ID}.original_mt.fa"

    if [[ "\${ROUND2_VCF}" == *.gz ]]; then
        cp "\${ROUND2_VCF}" "\${SAMPLE_ID}.round2.consensus_coord.input.vcf.gz"
    else
        bgzip -c "\${ROUND2_VCF}" > "\${SAMPLE_ID}.round2.consensus_coord.input.vcf.gz"
    fi
    tabix -p vcf -f "\${SAMPLE_ID}.round2.consensus_coord.input.vcf.gz"

    python3 - <<'PY_LIFTBACK'
import gzip
import re
from pathlib import Path
from collections import Counter

sample = "${meta.id}"
orig_fasta = Path(f"{sample}.original_mt.fa")
cons_sites_vcf = Path("${consensus_sites_vcf}")
round2_vcf = Path(f"{sample}.round2.consensus_coord.input.vcf.gz")
out_vcf = Path(f"{sample}.round2.original_coords.clean.final.split.vcf")
unresolved_tsv = Path(f"{sample}.round2.original_coords.unresolved.tsv")
skipped_overlap_tsv = Path(f"{sample}.round2.original_coords.skipped_overlap_consensus_sites.tsv")
summary_tsv = Path(f"{sample}.round2.original_coords.liftback.summary.tsv")


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def read_single_fasta(path):
    name = None
    seq_parts = []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    raise SystemExit(f"ERROR: {path} contains more than one contig; this lift-back expects one mt contig")
                name = line[1:].split()[0]
            else:
                seq_parts.append(line.strip())
    if name is None:
        raise SystemExit(f"ERROR: empty FASTA: {path}")
    return name, "".join(seq_parts).upper()


def parse_info(info):
    if info in (".", ""):
        return []
    return info.split(";")


def add_info(info, items):
    parts = [] if info in (".", "") else info.split(";")
    # Avoid duplicate keys if the WDL was rerun through this process.
    drop_keys = {x.split("=", 1)[0] for x in items}
    kept = [x for x in parts if x.split("=", 1)[0] not in drop_keys]
    return ";".join(kept + items) if kept or items else "."


def split_sample(fmt, sample_value):
    keys = fmt.split(":") if fmt not in (".", "") else []
    vals = sample_value.split(":") if sample_value not in (".", "") else []
    return keys, vals


def swap_two_value_array(value):
    if value in (".", ""):
        return value
    arr = value.split(",")
    if len(arr) == 2:
        return ",".join([arr[1], arr[0]])
    return value


def invert_af(value, ad_value=None):
    if value in (".", ""):
        return value
    # Prefer AD-derived AF when possible.
    if ad_value not in (None, ".", ""):
        try:
            ad = [float(x) for x in ad_value.split(",")]
            if len(ad) == 2 and sum(ad) > 0:
                return f"{ad[0] / sum(ad):.6g}"
        except Exception:
            pass
    try:
        vals = value.split(",")
        if len(vals) == 1:
            return f"{1.0 - float(vals[0]):.6g}"
    except Exception:
        pass
    return value


def flip_sample_fields(fmt, sample_value):
    keys, vals = split_sample(fmt, sample_value)
    if not keys or not vals:
        return sample_value
    if len(vals) < len(keys):
        vals = vals + ["."] * (len(keys) - len(vals))
    d = dict(zip(keys, vals))
    old_ad = d.get("AD")

    # Arrays ordered as REF,ALT.
    for key in ("AD", "F1R2", "F2R1"):
        if key in d:
            d[key] = swap_two_value_array(d[key])

    # Strand bias in Mutect-style VCF is REF_FWD,REF_REV,ALT_FWD,ALT_REV.
    if "SB" in d and d["SB"] not in (".", ""):
        sb = d["SB"].split(",")
        if len(sb) == 4:
            d["SB"] = ",".join([sb[2], sb[3], sb[0], sb[1]])

    if "AF" in d:
        d["AF"] = invert_af(d["AF"], old_ad)

    # Keep heteroplasmic genotype style stable after flipping.
    if "GT" in d and d["GT"] not in (".", "./."):
        sep = "/" if "/" in d["GT"] else "|" if "|" in d["GT"] else None
        if sep:
            alleles = d["GT"].split(sep)
            if set(alleles).issubset({"0", "1"}):
                d["GT"] = sep.join("1" if a == "0" else "0" for a in alleles)
                # For split biallelic mtDNA calls, 0/1 and 1/0 have the same meaning;
                # normalize to 0/1 to avoid confusing downstream tools.
                if set(d["GT"].replace("|", "/").split("/")) == {"0", "1"}:
                    d["GT"] = "0/1" if sep == "/" else "0|1"

    return ":".join(d.get(k, ".") for k in keys)


def parse_vcf_records(path):
    records = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\\n").split("\\t")
            if len(f) < 5:
                continue
            chrom, pos, _id, ref, alt = f[:5]
            if "," in alt:
                raise SystemExit(f"ERROR: consensus-sites VCF must be split/biallelic, found multi-ALT at {chrom}:{pos} {ref}>{alt}")
            if alt.startswith("<") or ref == "*" or alt == "*":
                raise SystemExit(f"ERROR: unsupported symbolic allele in consensus-sites VCF at {chrom}:{pos} {ref}>{alt}")
            records.append((int(pos), ref.upper(), alt.upper()))
    records.sort(key=lambda x: x[0])
    return records


orig_name, orig_seq = read_single_fasta(orig_fasta)
orig_len = len(orig_seq)
cons_events = parse_vcf_records(cons_sites_vcf)

# Build the same consensus sequence coordinate system as bcftools consensus -H A,
# while retaining a per-base map from consensus position back to original position.
cons_seq_parts = []
# map_pos is 1-based consensus-position indexed after adding a dummy at element 0.
map_pos = [None]
map_kind = [None]       # "orig" for original-backed base; "ins" for inserted base
map_anchor = [None]     # for inserted bases, original anchor position
cur0 = 0
cumulative_delta = 0
applied_consensus_events = 0
skipped_overlap_events = []

for pos1, ref, alt in cons_events:
    start0 = pos1 - 1
    end0 = start0 + len(ref)
    if start0 < cur0:
        # bcftools consensus skips records that overlap a previously applied variant.
        # To reconstruct the same coordinate system, skip the current overlapping event
        # rather than aborting the lift-back step.
        skipped_overlap_events.append((pos1, ref, alt, cur0))
        continue
    observed = orig_seq[start0:end0]
    if observed != ref:
        raise SystemExit(f"ERROR: REF mismatch in consensus-sites VCF at original {pos1}: VCF REF={ref}, FASTA={observed}")

    # unchanged original segment before this variant
    for i in range(cur0, start0):
        cons_seq_parts.append(orig_seq[i])
        map_pos.append(i + 1)
        map_kind.append("orig")
        map_anchor.append(i + 1)

    # ALT segment applied to consensus
    for i, base in enumerate(alt):
        cons_seq_parts.append(base)
        if i < len(ref):
            map_pos.append(pos1 + i)
            map_kind.append("orig")
            map_anchor.append(pos1 + i)
        else:
            # Inserted base after the first/anchor REF base.
            map_pos.append(None)
            map_kind.append("ins")
            map_anchor.append(pos1 + len(ref) - 1)

    cumulative_delta += len(alt) - len(ref)
    cur0 = end0
    applied_consensus_events += 1

for i in range(cur0, orig_len):
    cons_seq_parts.append(orig_seq[i])
    map_pos.append(i + 1)
    map_kind.append("orig")
    map_anchor.append(i + 1)

cons_seq = "".join(cons_seq_parts)
cons_len = len(cons_seq)

with skipped_overlap_tsv.open("w") as sk:
    sk.write("sample	original_pos	ref	alt	previous_applied_original_end0	reason\\n")
    for pos1, ref, alt, prev_end0 in skipped_overlap_events:
        sk.write(f"{sample}	{pos1}	{ref}	{alt}	{prev_end0}	skipped_to_match_bcftools_consensus_overlap_behavior\\n")

# Header lines.
header = []
with open_text(round2_vcf) as handle:
    for line in handle:
        if line.startswith("##"):
            if line.startswith("##contig=<ID="):
                # Replace contig length with the original reference length.
                m = re.match(r"##contig=<ID=([^,>]+).*", line.rstrip("\\n"))
                chrom_id = m.group(1) if m else orig_name
                header.append(f"##contig=<ID={chrom_id},length={orig_len},assembly={orig_fasta.name}>\\n")
            else:
                header.append(line)
        elif line.startswith("#CHROM"):
            header.append('##INFO=<ID=CONS_POS,Number=1,Type=Integer,Description="Position of this record on the sample-specific consensus reference before lift-back">\\n')
            header.append('##INFO=<ID=CONS_REF,Number=1,Type=String,Description="REF allele of this record on the sample-specific consensus reference before lift-back">\\n')
            header.append('##INFO=<ID=CONS_ALT,Number=1,Type=String,Description="ALT allele of this record on the sample-specific consensus reference before lift-back">\\n')
            header.append('##INFO=<ID=LIFTBACK_STATUS,Number=1,Type=String,Description="How the consensus-coordinate record was converted to original-reference coordinates">\\n')
            header.append(line)
            break

stats = Counter()
unresolved_rows = []

with out_vcf.open("w") as out, unresolved_tsv.open("w") as unres:
    out.writelines(header)
    unres.write("sample\\tcons_chrom\\tcons_pos\\tcons_ref\\tcons_alt\\treason\\tmapped_original_pos\\toriginal_ref_candidate\\n")

    with open_text(round2_vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            line = line.rstrip("\\n")
            if not line:
                continue
            f = line.split("\\t")
            if len(f) < 8:
                stats["unresolved_malformed"] += 1
                continue

            chrom, pos_s, vid, cref, calt, qual, flt, info = f[:8]
            pos = int(pos_s)
            cref = cref.upper()
            alts = [a.upper() for a in calt.split(",")]

            if len(alts) != 1:
                stats["unresolved_multiallelic"] += 1
                unres.write(f"{sample}\\t{chrom}\\t{pos}\\t{cref}\\t{calt}\\tmultiallelic_round2_record\\t.\\t.\\n")
                continue
            if pos < 1 or pos + len(cref) - 1 > cons_len:
                stats["unresolved_out_of_consensus_range"] += 1
                unres.write(f"{sample}\\t{chrom}\\t{pos}\\t{cref}\\t{calt}\\tout_of_consensus_range\\t.\\t.\\n")
                continue
            observed_cref = cons_seq[pos-1:pos-1+len(cref)]
            if observed_cref != cref:
                stats["unresolved_consensus_ref_mismatch"] += 1
                unres.write(f"{sample}\\t{chrom}\\t{pos}\\t{cref}\\t{calt}\\tconsensus_ref_mismatch_observed_{observed_cref}\\t.\\t.\\n")
                continue

            q_positions = range(pos, pos + len(cref))
            backed = [map_pos[q] for q in q_positions if map_pos[q] is not None]
            anchors = [map_anchor[q] for q in q_positions if map_anchor[q] is not None]
            if not backed and not anchors:
                stats["unresolved_no_original_anchor"] += 1
                unres.write(f"{sample}\\t{chrom}\\t{pos}\\t{cref}\\t{calt}\\tno_original_anchor\\t.\\t.\\n")
                continue

            # For inserted bases in the consensus, use the original anchor. For normal bases,
            # use the consecutive original span covered by the record.
            start_orig = min(backed) if backed else min(anchors)
            end_orig = max(backed) if backed else start_orig
            orig_ref_candidate = orig_seq[start_orig-1:end_orig]

            allele_alt = alts[0]
            status = None
            new_ref = None
            new_alt = None
            new_pos = start_orig
            flip = False

            if orig_ref_candidate == cref:
                status = "REF_MATCH_SHIFT_ONLY"
                new_ref = orig_ref_candidate
                new_alt = allele_alt
                flip = False
            elif orig_ref_candidate == allele_alt:
                status = "REF_ALT_FLIP"
                new_ref = orig_ref_candidate
                new_alt = cref
                flip = True
            else:
                # For pure insertion into the consensus, the original REF is often the anchor base,
                # which can match the ALT allele after flipping.
                if allele_alt == orig_ref_candidate:
                    status = "REF_ALT_FLIP"
                    new_ref = orig_ref_candidate
                    new_alt = cref
                    flip = True
                else:
                    stats["unresolved_original_ref_not_in_round2_alleles"] += 1
                    unres.write(f"{sample}\\t{chrom}\\t{pos}\\t{cref}\\t{calt}\\toriginal_ref_not_in_round2_alleles\\t{start_orig}\\t{orig_ref_candidate}\\n")
                    continue

            # Final REF sanity check against original FASTA.
            if orig_seq[new_pos-1:new_pos-1+len(new_ref)] != new_ref:
                stats["unresolved_final_ref_mismatch"] += 1
                unres.write(f"{sample}\\t{chrom}\\t{pos}\\t{cref}\\t{calt}\\tfinal_ref_mismatch\\t{new_pos}\\t{new_ref}\\n")
                continue

            f[0] = orig_name if orig_name else chrom
            f[1] = str(new_pos)
            f[3] = new_ref
            f[4] = new_alt
            f[7] = add_info(info, [f"CONS_POS={pos}", f"CONS_REF={cref}", f"CONS_ALT={allele_alt}", f"LIFTBACK_STATUS={status}"])

            if flip and len(f) >= 10:
                f[9] = flip_sample_fields(f[8], f[9])

            out.write("\\t".join(f) + "\\n")
            stats["resolved_total"] += 1
            stats[f"resolved_{status}"] += 1
            if len(new_ref) == 1 and len(new_alt) == 1:
                stats["resolved_snv"] += 1
            elif len(new_ref) < len(new_alt):
                stats["resolved_insertion_vs_original"] += 1
            elif len(new_ref) > len(new_alt):
                stats["resolved_deletion_vs_original"] += 1
            else:
                stats["resolved_complex_same_length"] += 1

with summary_tsv.open("w") as s:
    s.write("metric\\tvalue\\n")
    s.write(f"original_reference\\t{orig_fasta}\\n")
    s.write(f"original_length\\t{orig_len}\\n")
    s.write(f"consensus_length_reconstructed_from_sites\\t{cons_len}\\n")
    s.write(f"consensus_net_delta\\t{cons_len - orig_len}\\n")
    s.write(f"consensus_sites_raw_n\\t{len(cons_events)}\\n")
    s.write(f"consensus_sites_applied_n\\t{applied_consensus_events}\\n")
    s.write(f"consensus_sites_skipped_overlap_n\\t{len(skipped_overlap_events)}\\n")
    for k in sorted(stats):
        s.write(f"{k}\\t{stats[k]}\\n")

print(f"[INFO] Reconstructed consensus length: {cons_len}; original length: {orig_len}; net delta: {cons_len - orig_len}")
print(f"[INFO] Skipped overlapping consensus-site records: {len(skipped_overlap_events)}")
print(f"[INFO] Resolved records: {stats['resolved_total']}")
print(f"[INFO] Wrote unresolved records to: {unresolved_tsv}")
PY_LIFTBACK

    # Sort and left-normalize the lifted-back Round 2 records on the original mt reference.
    # This file contains only residual variants called relative to the consensus reference.
    bcftools sort \
        -Oz -o "\${SAMPLE_ID}.round2.original_coords.round2_only.unsorted.vcf.gz" \
        "\${SAMPLE_ID}.round2.original_coords.clean.final.split.vcf"

    bcftools index -t -f "\${SAMPLE_ID}.round2.original_coords.round2_only.unsorted.vcf.gz"

    bcftools norm \
        -f "\${SAMPLE_ID}.original_mt.fa" \
        -m-any \
        "\${SAMPLE_ID}.round2.original_coords.round2_only.unsorted.vcf.gz" \
        -Oz -o "\${SAMPLE_ID}.round2.original_coords.round2_only.clean.final.split.vcf.gz"

    tabix -p vcf -f "\${SAMPLE_ID}.round2.original_coords.round2_only.clean.final.split.vcf.gz"

    # Merge back the high-confidence Round 1 homoplasmies used to build the consensus.
    # Merge policy:
    #   1) Round 1 consensus-basis variants have priority.
    #   2) Keep all exact-unique Round 1 PASS homoplasmies, including Round1-vs-Round1
    #      overlapping variants; do not drop them simply because bcftools consensus may
    #      have skipped an overlapping site during consensus FASTA construction.
    #   3) If a Round 2 residual call is exactly the same as a Round 1 consensus-basis variant,
    #      keep only the Round 1 record.
    #   4) If a Round 2 residual call overlaps the genomic span of a Round 1 consensus-basis
    #      variant, drop the Round 2 record to avoid reverse/residual artifacts around
    #      consensus indels, e.g. CT>C plus C>CT at the same locus.
    #   5) Keep only non-overlapping Round 2 residual calls.
    # Normalize the Round 1 consensus-basis VCF first so exact duplicate and overlap checks
    # are performed in the same representation as the lifted Round 2 VCF.
    bcftools norm \
        -f "\${SAMPLE_ID}.original_mt.fa" \
        -m-any \
        "${consensus_sites_vcf}" \
        -Oz -o "\${SAMPLE_ID}.round1_consensus_basis.normalized.vcf.gz"

    tabix -p vcf -f "\${SAMPLE_ID}.round1_consensus_basis.normalized.vcf.gz"

    python3 - <<'PY_MERGE_CONSENSUS_SITES'
import gzip
import re
from pathlib import Path
from collections import Counter

sample = "${meta.id}"
round2_only_vcf = Path(f"{sample}.round2.original_coords.round2_only.clean.final.split.vcf.gz")
consensus_sites_vcf = Path(f"{sample}.round1_consensus_basis.normalized.vcf.gz")
merged_vcf = Path(f"{sample}.round2.original_coords.with_consensus_basis.unsorted.vcf")
conflict_tsv = Path(f"{sample}.round2.original_coords.round2_dropped_by_round1.tsv")
merge_summary_tsv = Path(f"{sample}.round2.original_coords.merge_policy.summary.tsv")


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def meta_key(line):
    line = line.rstrip("\\n")
    m = re.match(r"##(INFO|FORMAT|FILTER|ALT|contig)=<ID=([^,>]+).*", line)
    if m:
        return (m.group(1), m.group(2))
    if line.startswith("##") and "=" in line:
        return ("META", line.split("=", 1)[0])
    return ("LINE", line)


def parse_info(info):
    if info in (".", ""):
        return []
    return [x for x in info.split(";") if x]


def set_source(info, source):
    parts = [x for x in parse_info(info) if x.split("=", 1)[0] != "SOURCE"]
    parts.append(f"SOURCE={source}")
    return ";".join(parts) if parts else "."


def read_headers(path):
    meta = []
    chrom_line = None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                meta.append(line.rstrip("\\n"))
            elif line.startswith("#CHROM"):
                chrom_line = line.rstrip("\\n")
                break
    if chrom_line is None:
        raise SystemExit(f"ERROR: missing #CHROM header in {path}")
    return meta, chrom_line


def iter_records(path):
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            line = line.rstrip("\\n")
            if not line:
                continue
            f = line.split("\\t")
            if len(f) < 8:
                continue
            yield f


def record_key(f):
    return (f[0], int(f[1]), f[3], f[4])


def record_span(f):
    # Return the original-reference interval occupied by the VCF REF allele.
    chrom = f[0]
    start = int(f[1])
    end = start + len(f[3]) - 1
    return chrom, start, end


def overlaps(a, b):
    achr, astart, aend = a
    bchr, bstart, bend = b
    return achr == bchr and astart <= bend and bstart <= aend


def variant_string(f):
    return f"{f[0]}:{f[1]}:{f[3]}>{f[4]}"


r2_meta, r2_chrom = read_headers(round2_only_vcf)
cs_meta, cs_chrom = read_headers(consensus_sites_vcf)

# Build a combined header. Keep the Round 2 #CHROM/sample line as authoritative.
seen_keys = set()
combined_meta = []
for line in r2_meta:
    key = meta_key(line)
    seen_keys.add(key)
    combined_meta.append(line)

for line in cs_meta:
    # Do not import contig lines from the original consensus-sites VCF; the lifted Round 2
    # header already carries the original mt contig definition used downstream.
    if line.startswith("##contig=<ID="):
        continue
    key = meta_key(line)
    if key not in seen_keys:
        seen_keys.add(key)
        combined_meta.append(line)

if ("INFO", "SOURCE") not in seen_keys:
    combined_meta.append('##INFO=<ID=SOURCE,Number=1,Type=String,Description="Variant source after Round 2 lift-back: round2_consensus_call or round1_consensus_basis">')
if ("INFO", "ROUND2_DROP_REASON") not in seen_keys:
    combined_meta.append('##INFO=<ID=ROUND2_DROP_REASON,Number=1,Type=String,Description="Diagnostic reason a Round 2 residual call was removed during Round 1 priority merge">')

stats = Counter()

# Read Round 1 consensus-basis variants first because they have priority.
# Deduplicate exact duplicate Round 1 records after normalization.
consensus_records = []
consensus_keys = set()
consensus_spans = []
for f in iter_records(consensus_sites_vcf):
    key = record_key(f)
    if key in consensus_keys:
        stats["round1_exact_duplicate_records_skipped"] += 1
        continue

    cs_span = record_span(f)
    # Keep all Round 1 consensus-basis records after exact-deduplication.
    # Important: bcftools consensus may skip overlapping Round 1 records while building
    # the consensus FASTA, but for the final original-coordinate VCF we still want to
    # report every high-confidence Round 1 PASS homoplasmy used as consensus-basis input.
    # Therefore, do NOT drop Round 1 records merely because they overlap another Round 1
    # record.  Round 1 still has priority over Round 2: any Round 2 residual call that
    # exactly matches or overlaps any kept Round 1 record will be removed below.
    if any(overlaps(cs_span, kept_span) for kept_span, _ in consensus_spans):
        stats["round1_overlap_records_kept_in_final_vcf"] += 1

    f[7] = set_source(f[7], "round1_consensus_basis")
    consensus_keys.add(key)
    consensus_records.append(f)
    consensus_spans.append((cs_span, variant_string(f)))

# Keep only Round 2 calls that are not exact duplicates of, and do not overlap,
# any Round 1 consensus-basis variant.
round2_records = []
with conflict_tsv.open("w") as conflict:
    conflict.write("sample\\tround2_variant\\tdrop_reason\\toverlapping_round1_variant\\n")
    for f in iter_records(round2_only_vcf):
        key = record_key(f)
        r2_var = variant_string(f)

        if key in consensus_keys:
            stats["round2_exact_duplicate_dropped_round1_preferred"] += 1
            conflict.write(f"{sample}\\t{r2_var}\\texact_duplicate_round1_preferred\\t{r2_var}\\n")
            continue

        r2_span = record_span(f)
        hit = None
        for cs_span, cs_var in consensus_spans:
            if overlaps(r2_span, cs_span):
                hit = cs_var
                break

        if hit is not None:
            stats["round2_overlap_dropped_round1_preferred"] += 1
            conflict.write(f"{sample}\\t{r2_var}\\toverlap_round1_consensus_basis\\t{hit}\\n")
            continue

        f[7] = set_source(f[7], "round2_consensus_call")
        round2_records.append(f)

with merged_vcf.open("w") as out:
    for line in combined_meta:
        print(line, file=out)
    print(r2_chrom, file=out)

    # Put Round 1 first so if a downstream tool ever sees an equivalent representation,
    # the priority source is visible first. Final bcftools sort will order by coordinate.
    for f in consensus_records:
        print("\\t".join(f), file=out)
    for f in round2_records:
        print("\\t".join(f), file=out)

stats["round1_consensus_basis_records_kept"] = len(consensus_records)
stats["round2_consensus_call_records_kept_after_round1_priority"] = len(round2_records)
stats["merged_unsorted_records_n"] = len(consensus_records) + len(round2_records)

with merge_summary_tsv.open("w") as ms:
    ms.write("metric\\tvalue\\n")
    for k in sorted(stats):
        ms.write(f"{k}\\t{stats[k]}\\n")

print(f"[INFO] Round1 consensus-basis records kept: {len(consensus_records)}")
print(f"[INFO] Round2 residual records kept after Round1 priority filtering: {len(round2_records)}")
print(f"[INFO] Round2 exact duplicates dropped: {stats['round2_exact_duplicate_dropped_round1_preferred']}")
print(f"[INFO] Round2 overlapping records dropped: {stats['round2_overlap_dropped_round1_preferred']}")
print(f"[INFO] Merge conflict report: {conflict_tsv}")
PY_MERGE_CONSENSUS_SITES

    # Sort and normalize the merged final VCF. The output filename remains unchanged,
    # but now it contains both Round 2 residual calls and Round 1 consensus-basis homoplasmies.
    bcftools sort \
        -Oz -o "\${SAMPLE_ID}.round2.original_coords.with_consensus_basis.sorted.vcf.gz" \
        "\${SAMPLE_ID}.round2.original_coords.with_consensus_basis.unsorted.vcf"

    bcftools index -t -f "\${SAMPLE_ID}.round2.original_coords.with_consensus_basis.sorted.vcf.gz"

    bcftools norm \
        -f "\${SAMPLE_ID}.original_mt.fa" \
        -m-any \
        "\${SAMPLE_ID}.round2.original_coords.with_consensus_basis.sorted.vcf.gz" \
        -Oz -o "\${SAMPLE_ID}.round2.original_coords.clean.final.split.vcf.gz"

    tabix -p vcf -f "\${SAMPLE_ID}.round2.original_coords.clean.final.split.vcf.gz"

    {
      echo -e "round2_only_records_n	\$(bcftools view -H \${SAMPLE_ID}.round2.original_coords.round2_only.clean.final.split.vcf.gz | wc -l)"
      echo -e "final_with_consensus_basis_records_n	\$(bcftools view -H \${SAMPLE_ID}.round2.original_coords.clean.final.split.vcf.gz | wc -l)"
      echo -e "round1_consensus_basis_records_n	\$(bcftools view -H \${SAMPLE_ID}.round2.original_coords.clean.final.split.vcf.gz | grep -c 'SOURCE=round1_consensus_basis' || true)"
      echo -e "round2_consensus_call_records_n	\$(bcftools view -H \${SAMPLE_ID}.round2.original_coords.clean.final.split.vcf.gz | grep -c 'SOURCE=round2_consensus_call' || true)"
    } >> "\${SAMPLE_ID}.round2.original_coords.liftback.summary.tsv"

    if [[ -s "\${SAMPLE_ID}.round2.original_coords.merge_policy.summary.tsv" ]]; then
      tail -n +2 "\${SAMPLE_ID}.round2.original_coords.merge_policy.summary.tsv" >> "\${SAMPLE_ID}.round2.original_coords.liftback.summary.tsv"
    fi


    # -------------------------------------------------------------------------
    # Lift Round 2 per-base coverage from consensus chrM coordinates back to
    # original chrM coordinates. The WDL CoverageAtEveryBase task uses the
    # consensus self-reference, so its position column is not directly comparable
    # to original-reference coverage when consensus indels exist.
    # -------------------------------------------------------------------------
    python3 - <<'PY_FIND_COVERAGE' > round2_coverage_path.txt
from pathlib import Path
import sys

sample = "${meta.id}"
root = Path("${cromwell_dir}")

if not root.exists():
    print(f"ERROR: Missing Cromwell directory: {root}", file=sys.stderr)
    sys.exit(1)

preferred = f"{sample}.standard.chrM_assigned.bam.per_base_coverage.tsv"
candidates = []
for p in root.rglob("*.per_base_coverage.tsv"):
    if not p.is_file():
        continue
    name = p.name
    if name == preferred:
        candidates.append((0, p))
    elif "standard.chrM_assigned" in name and name.endswith(".per_base_coverage.tsv"):
        candidates.append((1, p))
    elif name.endswith(".per_base_coverage.tsv"):
        candidates.append((2, p))

if not candidates:
    print(f"ERROR: Cannot find WDL per-base coverage file under {root}", file=sys.stderr)
    for p in sorted(root.rglob("*CoverageAtEveryBase*")):
        print(f"[COVERAGE_PATH] {p}", file=sys.stderr)
    sys.exit(1)

candidates.sort(key=lambda x: (x[0], -x[1].stat().st_mtime))
if len(candidates) > 1:
    print("[WARN] Multiple per-base coverage candidates found. Selecting best match:", file=sys.stderr)
    for rank, p in candidates[:10]:
        print(f"  rank={rank} {p}", file=sys.stderr)

print(candidates[0][1])
PY_FIND_COVERAGE

    COVERAGE_FILE="\$(cat round2_coverage_path.txt)"
    echo "[INFO] Selected round2 consensus-coordinate coverage: \${COVERAGE_FILE}"

    cp "\${COVERAGE_FILE}" "\${SAMPLE_ID}.round2.consensus_coords.per_base_coverage.tsv"

    python3 - <<'PY_LIFT_COVERAGE'
import gzip
from pathlib import Path
from collections import Counter

sample = "${meta.id}"
orig_fasta = Path(f"{sample}.original_mt.fa")
cons_sites_vcf = Path("${consensus_sites_vcf}")
coverage_in = Path(f"{sample}.round2.consensus_coords.per_base_coverage.tsv")
coverage_out = Path(f"{sample}.round2.original_coords.per_base_coverage.tsv")
inserted_out = Path(f"{sample}.round2.original_coords.inserted_base_coverage.tsv")
deleted_out = Path(f"{sample}.round2.original_coords.deleted_positions_without_consensus_coverage.tsv")
summary_out = Path(f"{sample}.round2.original_coords.coverage_liftback.summary.tsv")


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def read_single_fasta(path):
    name = None
    seq_parts = []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    raise SystemExit(f"ERROR: {path} contains more than one contig")
                name = line[1:].split()[0]
            else:
                seq_parts.append(line.strip())
    if name is None:
        raise SystemExit(f"ERROR: empty FASTA: {path}")
    return name, "".join(seq_parts).upper()


def parse_vcf_records(path):
    records = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\\n").split("\\t")
            if len(f) < 5:
                continue
            chrom, pos, _id, ref, alt = f[:5]
            if "," in alt or alt.startswith("<") or ref == "*" or alt == "*":
                raise SystemExit(f"ERROR: unsupported consensus-site allele at {chrom}:{pos} {ref}>{alt}")
            records.append((int(pos), ref.upper(), alt.upper()))
    records.sort(key=lambda x: x[0])
    return records


orig_name, orig_seq = read_single_fasta(orig_fasta)
orig_len = len(orig_seq)
cons_events = parse_vcf_records(cons_sites_vcf)

# Reconstruct the original -> consensus coordinate map using the same behavior as
# bcftools consensus: overlapping consensus-site records are skipped.
cons_seq_parts = []
map_pos = [None]      # 1-based consensus coordinate -> original position, or None for inserted bases
map_kind = [None]     # orig or ins
map_anchor = [None]   # for inserted bases, original anchor position
cur0 = 0
stats = Counter()
skipped_overlap = 0

for pos1, ref, alt in cons_events:
    start0 = pos1 - 1
    end0 = start0 + len(ref)
    if start0 < cur0:
        skipped_overlap += 1
        continue
    observed = orig_seq[start0:end0]
    if observed != ref:
        raise SystemExit(f"ERROR: REF mismatch in consensus-sites VCF at original {pos1}: VCF REF={ref}, FASTA={observed}")

    for i in range(cur0, start0):
        cons_seq_parts.append(orig_seq[i])
        map_pos.append(i + 1)
        map_kind.append("orig")
        map_anchor.append(i + 1)

    for i, base in enumerate(alt):
        cons_seq_parts.append(base)
        if i < len(ref):
            map_pos.append(pos1 + i)
            map_kind.append("orig")
            map_anchor.append(pos1 + i)
        else:
            map_pos.append(None)
            map_kind.append("ins")
            map_anchor.append(pos1 + len(ref) - 1)

    cur0 = end0
    stats["consensus_events_applied_n"] += 1

for i in range(cur0, orig_len):
    cons_seq_parts.append(orig_seq[i])
    map_pos.append(i + 1)
    map_kind.append("orig")
    map_anchor.append(i + 1)

cons_len = len(cons_seq_parts)

if not coverage_in.exists() or coverage_in.stat().st_size == 0:
    raise SystemExit(f"ERROR: missing coverage input: {coverage_in}")

with coverage_in.open() as fh:
    first = fh.readline().rstrip("\\n")
    if not first:
        raise SystemExit(f"ERROR: empty coverage input: {coverage_in}")
    header = first.split("\\t")
    lower = [x.lower() for x in header]
    has_header = "pos" in lower or "position" in lower or "chrom" in lower or "chromosome" in lower

    if has_header:
        pos_col = lower.index("pos") if "pos" in lower else lower.index("position")
        chrom_col = lower.index("chrom") if "chrom" in lower else lower.index("chromosome") if "chromosome" in lower else 0
        data_iter = fh
        out_header = header + ["CONS_POS", "COVERAGE_LIFTBACK_STATUS"]
    else:
        # Fallback for simple samtools-depth-like files without a header.
        header = ["chrom", "pos", "depth"]
        chrom_col = 0
        pos_col = 1
        data_iter = [first] + list(fh)
        out_header = header + ["CONS_POS", "COVERAGE_LIFTBACK_STATUS"]

    lifted_rows = []
    inserted_rows = []
    observed_original_positions = set()

    for line in data_iter:
        line = line.rstrip("\\n")
        if not line:
            continue
        row = line.split("\\t")
        if len(row) <= pos_col:
            stats["coverage_rows_malformed"] += 1
            continue
        try:
            cons_pos = int(float(row[pos_col]))
        except ValueError:
            stats["coverage_rows_bad_position"] += 1
            continue

        if cons_pos < 1 or cons_pos >= len(map_pos):
            stats["coverage_rows_out_of_consensus_range"] += 1
            continue

        orig_pos = map_pos[cons_pos]
        if orig_pos is None:
            anchor = map_anchor[cons_pos]
            inserted_rows.append(row + [str(cons_pos), str(anchor), "inserted_base_in_consensus_no_original_coordinate"])
            stats["coverage_rows_inserted_consensus_base"] += 1
            continue

        row[chrom_col] = orig_name
        row[pos_col] = str(orig_pos)
        lifted_rows.append((orig_pos, cons_pos, row + [str(cons_pos), "orig_backed_base"]))
        observed_original_positions.add(orig_pos)
        stats["coverage_rows_lifted_to_original"] += 1

lifted_rows.sort(key=lambda x: (x[0], x[1]))

with coverage_out.open("w") as out:
    out.write("\\t".join(out_header) + "\\n")
    for _orig_pos, _cons_pos, row in lifted_rows:
        out.write("\\t".join(row) + "\\n")

with inserted_out.open("w") as out:
    out.write("\\t".join(header + ["CONS_POS", "ANCHOR_ORIGINAL_POS", "COVERAGE_LIFTBACK_STATUS"]) + "\\n")
    for row in inserted_rows:
        out.write("\\t".join(row) + "\\n")

all_original_positions = set(range(1, orig_len + 1))
deleted_positions = sorted(all_original_positions - observed_original_positions)
with deleted_out.open("w") as out:
    out.write("sample\\tchrom\\toriginal_pos\\toriginal_base\\treason\\n")
    for pos in deleted_positions:
        out.write(f"{sample}\\t{orig_name}\\t{pos}\\t{orig_seq[pos-1]}\\tno_consensus_base_maps_to_this_original_position_or_no_coverage_row\\n")

stats["original_length"] = orig_len
stats["consensus_length_reconstructed_from_sites"] = cons_len
stats["consensus_net_delta"] = cons_len - orig_len
stats["consensus_sites_raw_n"] = len(cons_events)
stats["consensus_sites_skipped_overlap_n"] = skipped_overlap
stats["deleted_or_missing_original_positions_n"] = len(deleted_positions)
stats["inserted_consensus_coverage_rows_n"] = len(inserted_rows)

with summary_out.open("w") as out:
    out.write("metric\\tvalue\\n")
    for k in sorted(stats):
        out.write(f"{k}\\t{stats[k]}\\n")

print(f"[INFO] Lifted coverage rows to original coordinates: {stats['coverage_rows_lifted_to_original']}")
print(f"[INFO] Inserted consensus-base coverage rows without original coordinate: {len(inserted_rows)}")
print(f"[INFO] Deleted/missing original positions in lifted coverage: {len(deleted_positions)}")
print(f"[INFO] Wrote lifted coverage: {coverage_out}")
PY_LIFT_COVERAGE

    echo "[INFO] Coverage lift-back summary:"
    cat "\${SAMPLE_ID}.round2.original_coords.coverage_liftback.summary.tsv"

    echo "[INFO] Lift-back summary:"
    cat "\${SAMPLE_ID}.round2.original_coords.liftback.summary.tsv"
    """
}
