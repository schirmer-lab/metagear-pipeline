process SUMMARIZE_MAG_CATALOG {
    tag "cohort"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.22/python_3.10.sif' :
        'docker.io/raphsoft/python_base:3.10-R4' }"

    // Cohort-level summary: one row per dRep cluster (cluster_id, winner,
    // n_members, n_member_samples, member_genomes, completeness, contamination,
    // gtdb_lineage). Single invocation per cohort.

    input:
    path cdb_csv
    path wdb_csv
    path genome_info_csv
    path gtdb_summary, stageAs: 'gtdb/*'

    output:
    path('mag_catalog.csv'), emit: catalog_summary
    path 'versions.yml',     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    summarize_mag_catalog.py \\
        --cdb ${cdb_csv} \\
        --wdb ${wdb_csv} \\
        --genome-info ${genome_info_csv} \\
        --gtdb-dir gtdb \\
        --out mag_catalog.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
