include { covermColumnSuffix } from './methods'

process COVERM_CONTIG {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::coverm==0.7.0--hb4818e0_2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coverm:0.7.0--hb4818e0_2' :
        'biocontainers/coverm:0.7.0--hb4818e0_2' }"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.abundance_*.tsv"), emit: abundance
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def requested = (task.ext.methods  ?: ['count']).collect { [ it, task.ext.args  ?: '' ] } +
                    (task.ext.methods2 ?: ['rpkm', 'tpm']).collect { [ it, task.ext.args2 ?: '' ] }
    def metric_calls = requested.collect { method, options ->
        "coverm contig --methods ${method} --bam-files ${bams} -t ${task.cpus} ${options} 1> ${prefix}.abundance_${method}.tsv 2> log_${method}.txt\n" +
        "    sed -i '1 s/ ${covermColumnSuffix(method)}//g' ${prefix}.abundance_${method}.tsv"
    }.join('\n\n    ')
    """
    ${metric_calls}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CoverM: \$(coverm --version | cut -d' ' -f2)
    END_VERSIONS
    """
}

process COVERM_CONTIG_MERGE {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.22/python_3.10.sif' :
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    tuple val(meta), path(tsv_files)

    output:
    tuple val(meta), path("*_merged.tsv"), emit: abundance_merged
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    coverm_merge.py ${tsv_files} -o ${prefix}_merged.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
