process PHOLD_COMPARE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::phold=1.2.5"
    // quay.io only — Galaxy-depot's SIF mirror is unreliable for this image
    // (~1.1 GB; routinely drops mid-pull). Singularity pulls from quay.io via
    // docker:// and converts to SIF locally with no Docker runtime needed.
    container 'quay.io/biocontainers/phold:1.2.5--pyhdfd78af_0'

    // Foldseek search of the 3Di-encoded queries against the PHOLD structural
    // reference DB. CPU-only — Foldseek's index scales with thread count, so
    // request a high cpu allocation here.
    //
    // Inputs:
    //   proteins_fasta — the AA FASTA used for prediction (PHOLD reads it
    //                    again to align Foldseek hits to the original sequence
    //                    for coverage/identity metrics)
    //   predict_dir    — the merged PHOLD_PREDICT output (phold_3di.fasta +
    //                    phold_aa.fasta over all shards). MERGE_PHOLD_PREDICTIONS
    //                    builds this upstream.
    //   phold_db       — directory from PHOLD_INSTALL (or params.phold_db).

    input:
    tuple val(meta), path(proteins_fasta), path(predict_dir)
    path phold_db

    output:
    tuple val(meta), path("compare_${meta.id}/phold_per_cds_predictions.tsv"), emit: per_cds_tsv
    tuple val(meta), path("compare_${meta.id}/phold_3di.fasta")              , emit: di_fasta, optional: true
    tuple val(meta), path("compare_${meta.id}/")                             , emit: compare_dir
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    phold proteins-compare \\
        -i ${proteins_fasta} \\
        --predictions_dir ${predict_dir} \\
        -o compare_${meta.id} \\
        -d ${phold_db} \\
        -t ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phold: \$(phold --version 2>&1 | sed 's/^phold, version //; s/^phold //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p compare_${meta.id}
    printf 'cds_id\\tphrog\\tfunction\\tproduct\\tbitscore\\tannotation_confidence\\n' \\
        > compare_${meta.id}/phold_per_cds_predictions.tsv
    touch compare_${meta.id}/phold_3di.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phold: 1.2.5
    END_VERSIONS
    """
}
