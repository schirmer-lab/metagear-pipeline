process MERGE_VIRAL_PHOLD {
    tag "cohort"
    label 'process_single'

    // Pure-stdlib python.
    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11':
        'biocontainers/python:3.11' }"

    // Joins PHOLD per-CDS annotations onto the viral catalog via the
    // viral_join_table (produced by FILTER_STRUCTURES_INPUTS).
    // Emits virus.proteins.phold.tsv keyed on viral_rep_id, with a
    // `relation` column recording whether each row's annotation came
    // direct, propagated through the protein cluster table, or is
    // missing (no annotation possible).

    input:
    tuple val(meta), path(phold_tsv), path(join_table)

    output:
    tuple val(meta), path("virus.proteins.phold.tsv"), emit: viral_phold
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    merge_viral_phold.py \\
        --phold-tsv ${phold_tsv} \\
        --join-table ${join_table} \\
        --out-tsv virus.proteins.phold.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'viral_rep_id\\trelation\\tphrog\\tfunction\\tproduct\\n' > virus.proteins.phold.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
