process PHAROKKA_PROTEINS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pharokka:1.7.3--pyhdfd78af_0':
        'biocontainers/pharokka:1.7.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(amino_fasta), path(contigs)
    path pharokka_db

    output:
    tuple val(meta), path("${prefix}_pharokka/*_proteins_full_merged_output.tsv")       , emit: merged_output
    tuple val(meta), path("${prefix}_pharokka/*_proteins_summary_output.tsv")           , emit: summary_output
    tuple val(meta), path("*-affi-contigs-for-dramv.tab")            , emit: affi_for_dramv
    path "versions.yml"                                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """

    INPUT=$amino_fasta
    if [[ $amino_fasta == *.gz ]]
    then
        gunzip -c $amino_fasta >| \$PWD/${prefix}_plain.fasta
        INPUT=\$PWD/${prefix}_plain.fasta
    fi

    pharokka_proteins.py \\
        --infile \${INPUT} \\
        --outdir ${prefix}_pharokka \\
        --d ${pharokka_db} \\
        --threads ${task.cpus} \\
        $args

    create_dram_input.py \\
        --pharokka ${prefix}_pharokka/*_proteins_summary_output.tsv \\
        --proteins-faa \${INPUT} \\
        --contigs-fasta ${contigs} \\
        --out-affi ${prefix}-affi-contigs-for-dramv.tab


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharokka: \$(pharokka.py --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}_pharokka
    touch ${prefix}_pharokka/${prefix}.gbk
    touch ${prefix}_pharokka/${prefix}.log
    touch ${prefix}_pharokka/${prefix}_cds_functions.tsv
    touch ${prefix}_pharokka/${prefix}_top_hits_card.tsv
    touch ${prefix}_pharokka/top_hits_vfdb.tsv
    touch ${prefix}_pharokka/${prefix}_top_hits_inphared
    touch ${prefix}_pharokka/${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharokka: \$(pharokka.py --version)
    END_VERSIONS
    """
}
