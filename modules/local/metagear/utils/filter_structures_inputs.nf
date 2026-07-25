process FILTER_STRUCTURES_INPUTS {
    tag "cohort"
    label 'process_low'

    // Pure-stdlib python — same container shape as CLASSIFY_GENES /
    // BATCH_BINS.
    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11':
        'biocontainers/python:3.11' }"

    // Selects the protein subset that will be fed through PHOLD's
    // ProstT5 + Foldseek pipeline, based on params.structures_scope.
    // Also emits the viral_join_table used by MERGE_VIRAL_PHOLD downstream
    // to populate virus.proteins.phold.tsv via cluster-membership lookups.
    //
    // The scope filter and viral top-up logic live in bin/filter_structures_inputs.py;
    // see that file for the precise rules per scope value.

    input:
    tuple val(meta), path(all_proteins_fasta), path(viral_proteins_fasta), path(clusters_tsv), path(pfam_tsv)

    output:
    tuple val(meta), path("input_subset.faa")     , emit: subset_fasta
    tuple val(meta), path("viral_join_table.tsv") , emit: viral_join_table
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    filter_structures_inputs.py \\
        --all-proteins-fasta ${all_proteins_fasta} \\
        --viral-proteins-fasta ${viral_proteins_fasta} \\
        --clusters-tsv ${clusters_tsv} \\
        --pfam-tsv ${pfam_tsv} \\
        --out-fasta input_subset.faa \\
        --out-join-tsv viral_join_table.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf '>stub_id\\nACGT\\n' > input_subset.faa
    printf 'viral_rep_id\\tphold_rep_id\\trelation\\nstub_id\\tstub_id\\tdirect\\n' \\
        > viral_join_table.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
