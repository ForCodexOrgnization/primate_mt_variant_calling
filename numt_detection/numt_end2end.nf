#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.numt_config = null

if (!params.numt_config) {
    error "Missing required parameter: --numt_config"
}

process RUN_NUMT_END2END {
    label 'numt_related'

    input:
    path numt_config

    output:
    path "numt_end2end.done"

    script:
    """
    bash ${projectDir}/run_numt_end2end.sh --config ${numt_config}
    touch numt_end2end.done
    """
}

workflow {
    RUN_NUMT_END2END(file(params.numt_config))
}
