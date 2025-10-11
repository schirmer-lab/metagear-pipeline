process DRAMV {
    tag "$meta.id"
    conda "bioconda::dram==1.5.0--pyhdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dram:1.5.0--pyhdfd78af_0' :
        'quay.io/biocontainers/dram:1.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(contigs), path(affi)
    path dram_db

    output:
    tuple val(meta), path("dramv-annotate/*.genes.faa"), emit: genes_faa
    tuple val(meta), path("dramv-annotate/*.genes.fna"), emit: genes_fna
    tuple val(meta), path("dramv-annotate/*.genes.gff"), emit: genes_gff
    tuple val(meta), path("dramv-annotate/*.annotations.tsv"), emit: annotations
    tuple val(meta), path("dramv-annotate/*.scaffolds.fna"), emit: scaffolds_fna
    tuple val(meta), path("dramv-distill/*.amg_summary.tsv"), emit: amg_summary
    tuple val(meta), path("dramv-distill/*.vMAG_stats.tsv"), emit: vmag_stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    INPUT=$contigs
    if [[ $contigs == *.gz ]]
    then
        gunzip -c $contigs >| \$PWD/${prefix}_plain.fasta
        INPUT=\$PWD/${prefix}_plain.fasta
    fi

    # Step 1: Ensure fully resolved absolute paths
    DB_PATH=\$(readlink -f ${dram_db})

    rm -rf dram_db_link
    mkdir -p dram_db_link
    ln -s "\$DB_PATH" dram_db_link/dram

    # Explicitly set environment variables with absolute paths
    export DRAM_DB_LOCATION=\$(readlink -f dram_db_link/dram)
    export DRAM_CONFIG_LOCATION=\${DRAM_DB_LOCATION}/CONFIG

    # Important fix: Clean any existing output dirs explicitly
    rm -rf dramv-annotate dramv-distill

    # Step 2: Run DRAM-v annotate
    DRAM-v.py annotate -i \${INPUT} -v ${affi} -o dramv-annotate ${args} --threads ${task.cpus}

    # step 2 summarize anntotations
    DRAM-v.py distill -i dramv-annotate/annotations.tsv -o dramv-distill

    # step 3 rename output files
    mv dramv-annotate/genes.faa dramv-annotate/${prefix}.genes.faa
    mv dramv-annotate/genes.fna dramv-annotate/${prefix}.genes.fna
    mv dramv-annotate/genes.gff dramv-annotate/${prefix}.genes.gff
    mv dramv-annotate/annotations.tsv dramv-annotate/${prefix}.annotations.tsv
    mv dramv-annotate/scaffolds.fna dramv-annotate/${prefix}.scaffolds.fna
    mv dramv-distill/amg_summary.tsv dramv-distill/${prefix}.amg_summary.tsv
    mv dramv-distill/vMAG_stats.tsv dramv-distill/${prefix}.vMAG_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        DRAM: 1.3
    END_VERSIONS
    """
}


