process FILTER_PRODIGAL {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.27/python_3.10.sif' :
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.filtered.fasta"), emit: filtered_fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_prodigal.py ${fasta} ${prefix}.filtered.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
