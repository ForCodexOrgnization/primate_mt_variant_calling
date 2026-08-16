#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.numt_config = null

if (!params.numt_config) {
    error "Missing required parameter: --numt_config"
}

def readShellConfigValue = { String configPath, String key ->
    def line = file(configPath).readLines().find { it ==~ /^\s*${java.util.regex.Pattern.quote(key)}\s*=.*/ }
    if (!line) {
        return null
    }
    def value = line.replaceFirst(/^\s*${java.util.regex.Pattern.quote(key)}\s*=\s*/, '').trim()
    if ((value.startsWith('"') && value.endsWith('"')) || (value.startsWith("'") && value.endsWith("'"))) {
        value = value.substring(1, value.length() - 1)
    }
    return value
}

def resolvedChrMRefDir = System.getenv('REF_DIR') ?: params.ref_dir
def resolvedWholeRefDir = System.getenv('GLOBAL_REF_DIR') ?: params.global_ref_dir
def resolvedNuclearOnlyRefDir = System.getenv('NUCLEAR_ONLY_REF_DIR') ?: params.nuclear_only_ref_dir

if (!resolvedChrMRefDir) {
    error "REF_DIR/params.ref_dir is not configured"
}
if (!resolvedWholeRefDir) {
    error "GLOBAL_REF_DIR/params.global_ref_dir is not configured"
}
if (!resolvedNuclearOnlyRefDir) {
    error "NUCLEAR_ONLY_REF_DIR/params.nuclear_only_ref_dir is not configured"
}

log.info "INFO: NUMT CHRM_REF_DIR=${resolvedChrMRefDir}"
log.info "INFO: NUMT WHOLE_REF_DIR=${resolvedWholeRefDir}"
log.info "INFO: NUMT NUCLEAR_ONLY_REF_DIR=${resolvedNuclearOnlyRefDir}"
log.info "INFO: NF_CONFIG_PATH=${System.getenv('NF_CONFIG_PATH') ?: System.getenv('NXF_CONFIG_FILE') ?: System.getenv('NF_CONFIG_FILE') ?: '<unknown>'}"

def configuredSample = readShellConfigValue(params.numt_config.toString(), 'SAMPLE') ?: readShellConfigValue(params.numt_config.toString(), 'SAMPLE_ID')
def configuredSpecies = readShellConfigValue(params.numt_config.toString(), 'SPECIES_NAME') ?: readShellConfigValue(params.numt_config.toString(), 'REF_NAME')
def configuredRef = readShellConfigValue(params.numt_config.toString(), 'REF_NAME') ?: configuredSpecies
def samplesTsv = readShellConfigValue(params.numt_config.toString(), 'SAMPLES_TSV') ?: readShellConfigValue(params.numt_config.toString(), 'SAMPLES_FILE')
if (!configuredSample && !samplesTsv) {
    error "Missing SAMPLE/SAMPLE_ID or SAMPLES_TSV/SAMPLES_FILE in ${params.numt_config}"
}

process RUN_NUMT_END2END {
    label 'numt_related'
    tag { sample_id }

    input:
    tuple val(sample_id), val(real_species), val(ref_species), val(config_mode), path(numt_config)

    output:
    path "${sample_id}.numt_end2end.done"

    script:
    def refSpecies = ref_species ?: real_species
    def launcherReferenceFingerprint = params.launcher_reference_fingerprint ?: 'unset'
    """
    set -euo pipefail
    echo "reference_fingerprint=${launcherReferenceFingerprint}" >/dev/null

    run_with_retry() {
      local cfg="\$1"
      local sample="\$2"
      local discovery_dir="\$3"
      local max_attempts="\${NUMT_MAX_ATTEMPTS:-2}"
      local attempt=1
      local status=0

      while (( attempt <= max_attempts )); do
        echo "[\$(date)] NUMT end-to-end attempt \${attempt}/\${max_attempts} for \${sample}"

        set +e
        bash ${projectDir}/run_numt_end2end.sh --config "\${cfg}"
        status=\$?
        set -e

        if [[ "\${status}" -eq 0 ]]; then
          return 0
        fi

        echo "[\$(date)] NUMT attempt \${attempt} failed with exit code \${status} for \${sample}" >&2

        if [[ ( "\${status}" -eq 137 || "\${status}" -eq 143 ) && "\${attempt}" -lt "\${max_attempts}" ]]; then
          echo "[\$(date)] Retrying \${sample} for possible transient/preemption failure." >&2
          rm -f "\${discovery_dir}/\${sample}.numt_discovery.done"
          rm -rf "\${discovery_dir}/intermediate" "\${discovery_dir}/tmp"
          attempt=\$((attempt + 1))
          sleep "\${NUMT_RETRY_SLEEP_SECONDS:-30}"
          continue
        fi

        return "\${status}"
      done

      return "\${status}"
    }

    # shellcheck disable=SC1090
    source ${numt_config}

    WHOLE_REF_DIR='${resolvedWholeRefDir}'
    NUCLEAR_ONLY_REF_DIR='${resolvedNuclearOnlyRefDir}'
    CHRM_REF_DIR='${resolvedChrMRefDir}'
    export WHOLE_REF_DIR NUCLEAR_ONLY_REF_DIR CHRM_REF_DIR

    SAMPLE_ID='${sample_id}'
    BESTHIT_DIR="\${BESTHIT_OUTDIR:-/nfs/roberts/pi/pi_njl27/lt692/primate_results_numt_besthit}"

    require_highconf_output() {
      local highconf_path="\${BESTHIT_DIR}/\${SAMPLE_ID}.highconf_numt.bed"
      [[ -e "\${highconf_path}" ]] || {
        echo "ERROR: Missing required file after NUMT end-to-end run: \${highconf_path}" >&2
        return 1
      }
    }

    if [[ '${config_mode}' == 'direct' ]]; then
      : "\${DISCOVERY_OUTDIR:?missing DISCOVERY_OUTDIR in config}"
      if run_with_retry ${numt_config} "\${SAMPLE_ID}" "\${DISCOVERY_OUTDIR}"; then
        require_highconf_output
        touch ${sample_id}.numt_end2end.done
        exit 0
      else
        status=\$?
        exit "\${status}"
      fi
    fi

    : "\${CRAM_ROOT_1:?missing CRAM_ROOT_1 in config}"
    : "\${CRAM_ROOT_2:?missing CRAM_ROOT_2 in config}"
    : "\${WHOLE_REF_DIR:?required config variable missing: WHOLE_REF_DIR}"
    : "\${NUCLEAR_ONLY_REF_DIR:?required config variable missing: NUCLEAR_ONLY_REF_DIR}"
    : "\${CHRM_REF_DIR:?required config variable missing: CHRM_REF_DIR}"
    : "\${DISCOVERY_OUTROOT:?missing DISCOVERY_OUTROOT in config}"
    : "\${BESTHIT_OUTDIR:?missing BESTHIT_OUTDIR in config}"

    resolve_cram_in_root() {
      local sid="\$1" root="\$2"
      local found
      found="\$(
        find "\${root}" -type f -name "\${sid}.cram" | head -n 1
      )"

      if [[ -z "\${found}" ]]; then
        found="\$(
          find "\${root}" -type f -name "\${sid}*.cram" | head -n 1
        )"
      fi

      printf '%s\n' "\${found}"
    }

    find_ref() {
      local dir="\$1" name="\$2" suffix
      shift 2
      for suffix in "\$@"; do
        if [[ -s "\${dir}/\${name}\${suffix}" ]]; then
          printf '%s\n' "\${dir}/\${name}\${suffix}"
          return 0
        fi
      done
      return 1
    }

    log_missing_ref() {
      local label="\$1" dir="\$2" name="\$3" suffix
      shift 3
      echo "ERROR: \${label} ref missing for \${name} under \${dir}" >&2
      echo "INFO: Candidate filenames checked:" >&2
      for suffix in "\$@"; do
        printf '%s\n' "\${dir}/\${name}\${suffix}" >&2
      done
      echo "INFO: Similar files found:" >&2
      find "\${dir}" -maxdepth 1 -iname "*\${name%%_*}*" -printf '%f\n' 2>/dev/null | head -n 20 >&2
    }

    REAL_SPECIES='${real_species}'
    REF_SPECIES='${refSpecies}'

    CRAM="\$(resolve_cram_in_root "\${SAMPLE_ID}" "\${CRAM_ROOT_1}")"
    [[ -n "\${CRAM}" ]] || CRAM="\$(resolve_cram_in_root "\${SAMPLE_ID}" "\${CRAM_ROOT_2}")"
    [[ -n "\${CRAM}" ]] || { echo "ERROR: CRAM not found for \${SAMPLE_ID}" >&2; exit 2; }

    if [[ -f "\${CRAM}.crai" ]]; then
      CRAI="\${CRAM}.crai"
    elif [[ -f "\${CRAM%.cram}.cram.crai" ]]; then
      CRAI="\${CRAM%.cram}.cram.crai"
    else
      echo "ERROR: CRAI not found for \${CRAM}" >&2
      exit 2
    fi

    WGS_REF="\$(find_ref "\${WHOLE_REF_DIR}" "\${REF_SPECIES}" ".fasta" ".fa" ".fna")" || { log_missing_ref "WGS" "\${WHOLE_REF_DIR}" "\${REF_SPECIES}" ".fasta" ".fa" ".fna"; exit 2; }
    NUCLEAR_REF="\$(find_ref "\${NUCLEAR_ONLY_REF_DIR}" "\${REF_SPECIES}" ".fasta" ".fa" ".fna" ".nuclear_only.fasta" ".nuclear_only.fa" ".nuclear_only.fna")" || { log_missing_ref "nuclear" "\${NUCLEAR_ONLY_REF_DIR}" "\${REF_SPECIES}" ".fasta" ".fa" ".fna" ".nuclear_only.fasta" ".nuclear_only.fa" ".nuclear_only.fna"; exit 2; }

    echo "INFO: sample=\${SAMPLE_ID}"
    echo "INFO: real_species=\${REAL_SPECIES}"
    echo "INFO: ref_species=\${REF_SPECIES}"
    echo "INFO: WGS_REF=\${WGS_REF}"
    echo "INFO: NUCLEAR_REF=\${NUCLEAR_REF}"

    samtools faidx "\${WGS_REF}" >/dev/null 2>&1 || true
    MT_CONTIG="\${MT_CONTIG:-chrM}"
    MT_LENGTH="\$(awk -v mt="\${MT_CONTIG}" '\$1==mt{print \$2; exit}' "\${WGS_REF}.fai")"
    [[ -n "\${MT_LENGTH}" ]] || { echo "ERROR: mt contig \${MT_CONTIG} not found in \${WGS_REF}.fai" >&2; exit 2; }

    SAMPLE_DISCOVERY_OUTDIR="\${DISCOVERY_OUTROOT}/\${SAMPLE_ID}"
    TMP_CFG="\$(mktemp "\${TMPDIR:-/tmp}/\${SAMPLE_ID}.numtcfg.XXXXXX")"
    cp ${numt_config} "\${TMP_CFG}"
    cat >> "\${TMP_CFG}" <<CFG
SAMPLE=\${SAMPLE_ID}
INPUT_BAM_CRAM=\${CRAM}
INPUT_INDEX=\${CRAI}
WGS_REF=\${WGS_REF}
NUCLEAR_REF=\${NUCLEAR_REF}
MT_CONTIG=\${MT_CONTIG}
MT_LENGTH=\${MT_LENGTH}
DISCOVERY_OUTDIR=\${SAMPLE_DISCOVERY_OUTDIR}
SAMPLES_TSV=\${SAMPLES_TSV:-}
WHOLE_REF_DIR=\${WHOLE_REF_DIR}
NUCLEAR_ONLY_REF_DIR=\${NUCLEAR_ONLY_REF_DIR}
CHRM_REF_DIR=\${CHRM_REF_DIR}
BESTHIT_OUTDIR=\${BESTHIT_OUTDIR}
CFG

    if run_with_retry "\${TMP_CFG}" "\${SAMPLE_ID}" "\${SAMPLE_DISCOVERY_OUTDIR}"; then
      rm -f "\${TMP_CFG}"
      require_highconf_output
      touch ${sample_id}.numt_end2end.done
      exit 0
    else
      status=\$?
      rm -f "\${TMP_CFG}"
      exit "\${status}"
    fi
    """
}

workflow {
    if (configuredSample) {
        ch_samples = Channel.of(tuple(configuredSample, configuredSpecies ?: configuredSample, configuredRef ?: configuredSpecies ?: configuredSample, 'direct', file(params.numt_config)))
    } else {
        ch_samples = Channel
            .fromPath(samplesTsv)
            .splitCsv(header: false, sep: '\t', strip: true)
            .filter { row -> row.size() >= 2 && row[0] && !row[0].startsWith('#') && !['sample', 'sample_id'].contains(row[0].toLowerCase()) }
            .map { row -> tuple(row[0], row[1], row.size() >= 3 && row[2] ? row[2] : row[1], 'batch', file(params.numt_config)) }
    }

    RUN_NUMT_END2END(ch_samples)
}
