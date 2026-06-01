process SKANI_TRIANGLE {
    tag "cohort"
    label 'process_high'

    // skani ships inside the dRep biocontainer (dRep v3.6+ uses skani as its
    // default secondary clustering backend). Reusing the dRep image keeps the
    // batching step container-free of new dependencies.
    conda "bioconda::drep=3.6.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/drep:3.6.2--pyhdfd78af_0':
        'quay.io/biocontainers/drep:3.6.2--pyhdfd78af_0' }"

    // Builds the cohort-wide pairwise ANI matrix used to define dRep batches.
    // skani triangle is the bulk all-vs-all comparison mode — much faster than
    // calling skani dist per pair. The `--sparse` form only emits pairs above
    // a coarse ANI floor (we pass --min-af 0 + --ci so we keep more edges
    // than the default cutoff and let BATCH_BINS apply the actual threshold
    // downstream).
    //
    // Output (TSV): Ref_file<TAB>Query_file<TAB>ANI<TAB>Align_fraction_ref<TAB>Align_fraction_query
    //               (header line + one row per kept pair)

    input:
    tuple val(meta), path(fastas, stageAs: 'input_fastas/*')

    output:
    tuple val(meta), path("cohort.skani.ani.tsv"), emit: ani
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    # List bin FASTAs by absolute path; skani triangle accepts -l <file>.
    find -L input_fastas/ -type f > bin_list.txt

    skani triangle \\
        -t ${task.cpus} \\
        --sparse \\
        -l bin_list.txt \\
        -o cohort.skani.ani.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version 2>&1 | sed -E 's/^skani //')
    END_VERSIONS
    """

    stub:
    """
    printf 'Ref_file\\tQuery_file\\tANI\\tAlign_fraction_ref\\tAlign_fraction_query\\n' > cohort.skani.ani.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: 0.2.2
    END_VERSIONS
    """
}
