process POST_MSPMINER {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.10 conda-forge::pandas conda-forge::click conda-forge::pyfaidx conda-forge::scikit-learn"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.27/python_3.10.sif' :
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    tuple val(meta), path(gene_catalog), path(all_msps)

    output:
    tuple val(meta), path("pangenome_sequences/"), emit: pangenome_dir
    tuple val(meta), path("pangenome_sequences/*.pangenome.fasta"), emit: pangenome_files
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Create output directory for pangenome sequences
    mkdir -p pangenome_sequences

    # Run the post-MSP mining utility to extract pangenome sequences
    postminer_utils.py helper get-msp-pangenome \\
        --gene-catalog-fp ${gene_catalog} \\
        --all-msps-fp ${all_msps} \\
        --msp-pangenome-dir pangenome_sequences/ \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pyfaidx: \$(python -c "import pyfaidx; print(pyfaidx.__version__)")
    END_VERSIONS
    """
}
