process BATCH_BINS {
    tag "cohort"
    label 'process_low'

    // Pure-stdlib python — same container shape as CLASSIFY_GENES /
    // MERGE_CONTIG_CLASSIFICATION.
    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11':
        'biocontainers/python:3.11' }"

    // Partitions cohort bins into dRep batches via connected components on
    // the sparse skani ANI TSV. See bin/batch_bins.py for the algorithm and
    // its safety argument (single-linkage at 0.90 strictly contains dRep's
    // 0.95 secondary clustering threshold).
    //
    // For each component the script produces a self-contained batch dir:
    //   batches/batch_<NNN>/
    //     bins/<bin>.fa              # hardlinks (or symlinks) to the bins
    //     drep_work_seed/            # mirrors STAGE_DREP_WORK's layout so
    //       genomeInfo.csv           #   the per-batch dRep call reuses the
    //                                #   same `--genomeInfo` ext.args path
    //                                #   (`drep_work/drep_work_seed/genomeInfo.csv`).
    //
    // The downstream subworkflow flattens the batch_<NNN>/ glob and fans
    // each entry out into its own DREP_DEREPLICATE invocation.

    input:
    tuple val(meta), path(ani_tsv), path(genome_info)
    path(bin_fastas, stageAs: 'input_bins/*')

    output:
    tuple val(meta), path("batches/batch_*", type: 'dir'), emit: batches
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    batch_bins.py \\
        --ani-tsv ${ani_tsv} \\
        --genome-info ${genome_info} \\
        --bins input_bins/*.fa \\
        --out-root batches \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p batches/batch_001/bins batches/batch_001/drep_work_seed
    printf 'genome,completeness,contamination\\nstub.fa,90,5\\n' > batches/batch_001/drep_work_seed/genomeInfo.csv
    touch batches/batch_001/bins/stub.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
