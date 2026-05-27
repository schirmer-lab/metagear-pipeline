process ENRICH_PER_CONTIG_TSV {
    tag "$meta.id"
    label 'process_single'

    // Pure-python (stdlib) join; reuse the python_base image used elsewhere
    // in the repo so we don't introduce a new container.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.22/python_3.10.sif' :
        'docker.io/raphsoft/python_base:3.10-R4' }"

    // Per-sample enrichment of v1's contigs.tsv with cohort cluster IDs and
    // GTDB-Tk lineage. Inputs are mostly cohort-global; only contigs_tsv and
    // contig_to_bin are per-sample.

    input:
    tuple val(meta), path(contigs_tsv), path(contig_to_bin)
    path cdb_csv
    path wdb_csv
    path gtdb_summary, stageAs: 'gtdb/*'

    output:
    tuple val(meta), path("${meta.id}.contigs.enriched.tsv"), emit: enriched
    path 'versions.yml',                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    enrich_per_contig_tsv.py \\
        --contigs ${contigs_tsv} \\
        --contig-to-bin ${contig_to_bin} \\
        --cdb ${cdb_csv} \\
        --wdb ${wdb_csv} \\
        --gtdb-dir gtdb \\
        --sample ${meta.id} \\
        --out ${meta.id}.contigs.enriched.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
