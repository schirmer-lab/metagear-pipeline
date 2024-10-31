process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.22/python_3.10.sif' :
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    path samplesheet
    val input_type

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in schirmer-lab/metagear-pipeline/bin/
    """
    input_validator.py \\
        --input $samplesheet \\
        --validation_type $input_type \\
        --output samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

process RENAME_FILES {

    label 'process_single'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.22/python_3.10.sif' :
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    tuple val(meta), val(position), path(files)

    output:
    tuple val(meta), val(position), path("*__*"), emit: renamed_files
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    tag = meta.tag ?: meta.id

    """
    for FILE in $files; do
        ln -s \$FILE ${tag}__\$(basename \$FILE)
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
