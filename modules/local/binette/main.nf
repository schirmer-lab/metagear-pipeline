process BINETTE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::binette=1.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/binette:1.2.1--pyh106432d_1':
        'quay.io/biocontainers/binette:1.2.1--pyh106432d_1' }"

    input:
    // Binette 1.2.1's CLI requires --bin_dirs XOR --contig2bin_tables (not both).
    // Both upstream binners (SemiBin2 + MetaBAT2) emit gzipped bin FASTAs in
    // directories, so we stage each into its own subdir and pass them as a
    // space-separated list to --bin_dirs.
    tuple val(meta), path(contigs), path(semibin_bins, stageAs: 'semibin_bins/*'), path(metabat_bins, stageAs: 'metabat_bins/*')
    path checkm2_db

    output:
    tuple val(meta), path("${prefix}/final_bins/*.fa")                       , emit: bins      , optional: true
    tuple val(meta), path("${prefix}/final_bins")                            , emit: bins_dir  , optional: true
    tuple val(meta), path("${prefix}/final_bins_quality_reports.tsv")        , emit: quality
    tuple val(meta), path("${prefix}")                                        , emit: outdir
    path "versions.yml"                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    binette \\
        --contigs ${contigs} \\
        --bin_dirs semibin_bins metabat_bins \\
        --checkm2_db ${checkm2_db} \\
        --outdir ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: \$(binette --version 2>&1 | sed 's/.*[Bb]inette //; s/[,\\s].*//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/final_bins
    touch ${prefix}/final_bins/bin_1.fa
    touch ${prefix}/final_bins_quality_reports.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: 1.2.1
    END_VERSIONS
    """
}
