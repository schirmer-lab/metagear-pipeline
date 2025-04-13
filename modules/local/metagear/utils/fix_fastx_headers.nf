process FIX_FASTX_HEADERS {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_val_*.re.fq.gz"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip -c ${reads[0]} > ${prefix}_val_1.fq
    gunzip -c ${reads[1]} > ${prefix}_val_2.fq

    fix_header.py ${prefix}_val_1.fq ${prefix}_val_1.re.fq /1
    fix_header.py ${prefix}_val_2.fq ${prefix}_val_2.re.fq /2

    gzip ${prefix}_val_1.re.fq
    gzip ${prefix}_val_2.re.fq

    rm ${prefix}_val_1.fq ${prefix}_val_2.fq
    """
}
