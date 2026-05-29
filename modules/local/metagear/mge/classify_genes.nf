process CLASSIFY_GENES {
    tag "$meta.id"
    label 'process_low'

    // Pure-python helper (bin/classify_genes.py); container only needs
    // Python 3 + the stdlib. Same container as MERGE_CONTIG_CLASSIFICATION.
    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11':
        'biocontainers/python:3.11' }"

    input:
    // per_contig_tsvs — list of per-sample <sample>.contigs.tsv from
    //                   MERGE_CONTIG_CLASSIFICATION (collected across samples)
    // clusters_tsv    — full gene clusters TSV (rep<TAB>member) from
    //                   gene_analysis (catalogs/genes/all.genes.clusters.tsv).
    //
    // Emits the canonical filename `all.genes.clusters.classified.tsv`. The
    // publishDir in conf/metagear/integrated_classification.config rewrites
    // it to `all.genes.clusters.classified.refined.tsv` at publish time to
    // match the .raw / .draft / .refined naming series used in classification/.
    tuple val(meta), path(per_contig_tsvs, stageAs: 'per_contig/*'), path(clusters_tsv)

    output:
    tuple val(meta), path("all.genes.clusters.classified.tsv"), emit: classified
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    classify_genes.py \\
        --per-contig-tsvs per_contig/*.contigs.tsv \\
        --clusters-tsv ${clusters_tsv} \\
        --output all.genes.clusters.classified.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'rep_id\\tclasses\\tnum_members\\tclass_counts\\trep_class\\tmulti_class\\n' > all.genes.clusters.classified.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
