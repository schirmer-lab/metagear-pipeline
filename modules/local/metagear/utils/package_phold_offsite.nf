process PACKAGE_PHOLD_OFFSITE {
    tag "cohort"
    label 'process_low'

    // Pure-stdlib python — same container shape as the other metagear utils.
    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11':
        'biocontainers/python:3.11' }"

    // Assembles the offsite-predict bundle for the structures workflow:
    // gathers per-shard FASTAs from SEQKIT_SPLIT2 + the viral_join_table.tsv
    // from FILTER_STRUCTURES_INPUTS into a self-contained directory the user
    // rsyncs to a GPU server, runs the bundled script, and rsyncs results back.
    //
    // The templates live under assets/structures/offsite/ and get rendered
    // with cohort-specific values (shard count, suggested walltime, ETA
    // estimates per GPU class). Output is one directory the structures.nf
    // subworkflow then publishes to ${outdir}/structures/offsite_predict/.

    input:
    tuple val(meta), path(shards, stageAs: 'shards_in/*'), path(viral_join_table)
    path readme_template
    path script_template
    path sbatch_template

    output:
    tuple val(meta), path("offsite_predict/"), emit: bundle
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def scope = task.ext.structures_scope ?: 'unknown'
    """
    package_phold_offsite.py \\
        --shards-dir shards_in \\
        --viral-join-table ${viral_join_table} \\
        --readme-template ${readme_template} \\
        --script-template ${script_template} \\
        --sbatch-template ${sbatch_template} \\
        --out-dir offsite_predict \\
        --structures-scope '${scope}' \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p offsite_predict/shards
    touch offsite_predict/README.md
    touch offsite_predict/run_phold_predict.sh
    touch offsite_predict/submit_phold_predict.sbatch
    touch offsite_predict/viral_join_table.tsv
    echo '{"n_shards": 0}' > offsite_predict/manifest.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
