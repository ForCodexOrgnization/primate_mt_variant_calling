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
    """
    set -euo pipefail

    run_with_retry() {
      local cfg="\$1"
      local sample="\$2"
      local discovery_dir="\$3"
      local max_attempts="\${NUMT_MAX_ATTEMPTS:-2}"
      local attempt=1
      local status=0

      while (( attempt <= max_attempts )); do
        echo "[\$(date)] NUMT end-to-end attempt \${attempt}/\${max_attempts} for \${sample}"
        if bash ${projectDir}/run_numt_end2end.sh --config "\${cfg}"; then
          status=0
          break
        fi

        status=\$?
        if [[ "\${status}" -eq 1 && "\${attempt}" -lt "\${max_attempts}" ]]; then
          echo "[\$(date)] NUMT attempt \${attempt} failed with exit code 1; retrying \${sample} once for possible transient filesystem failure." >&2
          rm -f "\${discovery_dir}/\${sample}.numt_discovery.done"
          rm -rf "\${discovery_dir}/intermediate" "\${discovery_dir}/tmp"
          attempt=\$((attempt + 1))
          sleep "\${NUMT_RETRY_SLEEP_SECONDS:-30}"
          continue
        fi

        break
      done

      return "\${status}"
    }

    # shellcheck disable=SC1090
    source ${numt_config}

    SAMPLE_ID='${sample_id}'

    if [[ '${config_mode}' == 'direct' ]]; then
      : "\${DISCOVERY_OUTDIR:?missing DISCOVERY_OUTDIR in config}"
      if run_with_retry ${numt_config} "\${SAMPLE_ID}" "\${DISCOVERY_OUTDIR}"; then
        touch ${sample_id}.numt_end2end.done
        exit 0
      else
        status=\$?
        exit "\${status}"
      fi
    fi

    : "\${CRAM_ROOT_1:?missing CRAM_ROOT_1 in config}"
    : "\${CRAM_ROOT_2:?missing CRAM_ROOT_2 in config}"
    : "\${WHOLE_REF_DIR:?missing WHOLE_REF_DIR in config}"
    : "\${NUCLEAR_ONLY_REF_DIR:?missing NUCLEAR_ONLY_REF_DIR in config}"
    : "\${CHRM_REF_DIR:?missing CHRM_REF_DIR in config}"
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

    REAL_SPECIES='${real_species}'
    REF_SPECIES='${refSpecies}'

    CRAM="\$(resolve_cram_in_root "\${SAMPLE_ID}" "\${CRAM_ROOT_1}")"
    [[ -n "\${CRAM}" ]] || CRAM="\$(resolve_cram_in_root "\${SAMPLE_ID}" "\${CRAM_ROOT_2}")"
    [[ -n "\${CRAM}" ]] || { echo "ERROR: CRAM not found for \${SAMPLE_ID}" >&2; exit 1; }

    if [[ -f "\${CRAM}.crai" ]]; then
      CRAI="\${CRAM}.crai"
    elif [[ -f "\${CRAM%.cram}.cram.crai" ]]; then
      CRAI="\${CRAM%.cram}.cram.crai"
    else
      echo "ERROR: CRAI not found for \${CRAM}" >&2
      exit 1
    fi

    WGS_REF="\$(find_ref "\${WHOLE_REF_DIR}" "\${REF_SPECIES}" ".fasta" ".fa" ".fna")" || { echo "ERROR: WGS ref missing for \${REF_SPECIES} under \${WHOLE_REF_DIR}" >&2; exit 1; }
    NUCLEAR_REF="\$(find_ref "\${NUCLEAR_ONLY_REF_DIR}" "\${REF_SPECIES}" ".fasta" ".fa" ".fna" ".nuclear_only.fasta" ".nuclear_only.fa" ".nuclear_only.fna")" || { echo "ERROR: nuclear ref missing for \${REF_SPECIES} under \${NUCLEAR_ONLY_REF_DIR}" >&2; exit 1; }

    samtools faidx "\${WGS_REF}" >/dev/null 2>&1 || true
    MT_CONTIG="\${MT_CONTIG:-chrM}"
    MT_LENGTH="\$(awk -v mt="\${MT_CONTIG}" '\$1==mt{print \$2; exit}' "\${WGS_REF}.fai")"
    [[ -n "\${MT_LENGTH}" ]] || { echo "ERROR: mt contig \${MT_CONTIG} not found in \${WGS_REF}.fai" >&2; exit 1; }

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
CHRM_REF_DIR=\${CHRM_REF_DIR}
BESTHIT_OUTDIR=\${BESTHIT_OUTDIR}
CFG

    if run_with_retry "\${TMP_CFG}" "\${SAMPLE_ID}" "\${SAMPLE_DISCOVERY_OUTDIR}"; then
      status=0
    else
      status=\$?
    fi
    rm -f "\${TMP_CFG}"
    [[ "\${status}" -eq 0 ]] || exit "\${status}"
    touch ${sample_id}.numt_end2end.done
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
