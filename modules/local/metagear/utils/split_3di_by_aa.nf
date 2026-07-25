process SPLIT_3DI_BY_AA {
    tag "$meta.id"
    label 'process_low'

    // Pure-stdlib python — same container shape as the other metagear utils.
    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11':
        'biocontainers/python:3.11' }"

    // Per-chunk scatter helper for PHOLD_COMPARE. SEQKIT_SPLIT2_COMPARE divides
    // the cohort AA FASTA into M chunks; this process takes one chunk + the
    // cohort-wide MERGE_PHOLD_PREDICTIONS output and emits a per-chunk
    // predictions_dir containing only that chunk's matching phold_aa.fasta
    // and phold_3di.fasta entries.
    //
    // The cohort merged_predict_dir is shared across all M scatter tasks
    // (Nextflow stages it per task; the symlink is cheap and the filter is
    // single-pass over each file).

    input:
    tuple val(meta), path(chunk_aa), path(merged_predict_dir)

    output:
    tuple val(meta), path("chunk_predict_dir/"), emit: chunk_predict_dir
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    split_3di_by_aa.py \\
        --chunk-aa ${chunk_aa} \\
        --merged-predict-dir ${merged_predict_dir} \\
        --out-dir chunk_predict_dir

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python3 --version 2>&1 | sed 's/Python //')
END_VERSIONS
    """

    stub:
    """
    mkdir -p chunk_predict_dir
    touch chunk_predict_dir/phold_aa.fasta
    touch chunk_predict_dir/phold_3di.fasta

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: 3.11.0
END_VERSIONS
    """
}
