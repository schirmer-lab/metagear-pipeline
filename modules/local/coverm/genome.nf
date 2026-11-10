include { covermColumnSuffix } from './methods'

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
    tuple val(meta), path("*.abundance_*.tsv"), emit: abundance
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def requested = (task.ext.methods  ?: ['count']).collect { [ it, task.ext.args  ?: '' ] } +
                    (task.ext.methods2 ?: ['rpkm', 'tpm']).collect { [ it, task.ext.args2 ?: '' ] }
    def metric_calls = requested.collect { method, options ->
        "coverm genome --methods ${method} --bam-files ${bams} --genome-definition ${genome_definition} -t ${task.cpus} ${options} 1> ${prefix}.abundance_${method}.tsv 2> log_${method}.txt\n" +
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
