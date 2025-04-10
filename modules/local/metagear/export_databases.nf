process EXPORT_DATABASES {

    tag "$meta.id"

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.22/python_3.10.sif' :
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    tuple val(meta), path(database), path(destination)

    output:
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    FULL_SRC=\$(readlink ${database})
    FULL_DESTINATION=\$(readlink ${destination})
    mkdir -p \$FULL_DESTINATION

    # Clean destination and copy database files
    rm -rf \$FULL_DESTINATION/*

    (cd \$FULL_SRC && tar cf - .) | (cd \$FULL_DESTINATION && tar xf -)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}