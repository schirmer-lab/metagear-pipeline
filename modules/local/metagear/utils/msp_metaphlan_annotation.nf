process MSP_METAPHLAN_ANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.10 conda-forge::pandas conda-forge::numpy conda-forge::scikit-learn"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/25.08.11/python_3.10.sif':
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    tuple val(meta), path(msp_profile), path(metaphlan_profile)
    val metaphlan_version

    output:
    tuple val(meta), path("msp_metaphlan*_LM.full.txt"), emit: full_results
    tuple val(meta), path("msp_metaphlan*_LM.bestR2.txt"), emit: best_results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def version = metaphlan_version ?: 'v3'
    """
    # Run MSP-MetaPhlAn taxonomic annotation using linear regression
    msp_metaphlan_LM.py \\
        --msp-profile ${msp_profile} \\
        --metaphlan-profile ${metaphlan_profile} \\
        --output-dir . \\
        --metaphlan-version ${version} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        scikit-learn: \$(python -c "import sklearn; print(sklearn.__version__)")
    END_VERSIONS
    """
}
