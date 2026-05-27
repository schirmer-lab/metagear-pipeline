process COVERM_GENOME {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::coverm==0.7.0--hb4818e0_2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coverm:0.7.0--hb4818e0_2' :
        'biocontainers/coverm:0.7.0--hb4818e0_2' }"

    // Parallel to COVERM_CONTIG, but invokes `coverm genome` with a
    // --genome-definition TSV that maps contigs back to their source MAG.
    // Output naming matches COVERM_CONTIG so the downstream merge module
    // (COVERM_CONTIG_MERGE) and its python helper (coverm_merge.py) work
    // unchanged for both modes.

    input:
    tuple val(meta), path(bams), path(genome_definition)

    output:
    tuple val(meta), path("*.abundance_count.tsv"), emit: abundance_count
    tuple val(meta), path("*.abundance_trimmed_mean.tsv"), emit: abundance_trimmed_mean
    tuple val(meta), path("*.abundance_rpkm.tsv"), emit: abundance_rpkm
    tuple val(meta), path("*.abundance_tpm.tsv"), emit: abundance_tpm
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ''
    def args2  = task.ext.args2  ?: ''
    """
    coverm genome --methods count --bam-files $bams --genome-definition ${genome_definition} -t $task.cpus $args 1> ${prefix}.abundance_count.tsv 2> log_count.txt
    sed -i '1 s/ Read Count//g' ${prefix}.abundance_count.tsv

    coverm genome --methods trimmed_mean --bam-files $bams --genome-definition ${genome_definition} -t $task.cpus $args2 1> ${prefix}.abundance_trimmed_mean.tsv 2> log_trimmed_mean.txt
    sed -i '1 s/ Trimmed Mean//g' ${prefix}.abundance_trimmed_mean.tsv

    coverm genome --methods rpkm --bam-files $bams --genome-definition ${genome_definition} -t $task.cpus $args2 1> ${prefix}.abundance_rpkm.tsv 2> log_rpkm.txt
    sed -i '1 s/ RPKM//g' ${prefix}.abundance_rpkm.tsv

    coverm genome --methods tpm --bam-files $bams --genome-definition ${genome_definition} -t $task.cpus $args2 1> ${prefix}.abundance_tpm.tsv 2> log_tpm.txt
    sed -i '1 s/ TPM//g' ${prefix}.abundance_tpm.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CoverM: \$(coverm --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
