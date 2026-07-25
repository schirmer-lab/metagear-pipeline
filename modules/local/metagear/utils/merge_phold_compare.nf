process MERGE_PHOLD_COMPARE {
    tag "cohort"
    label 'process_single'

    // coreutils container is sufficient (cat + awk only).
    conda "conda-forge::coreutils"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5':
        'biocontainers/coreutils:9.5' }"

    // Gather step for the PHOLD_COMPARE scatter. Each scatter task emits one
    // phold_per_cds_predictions.tsv (rows keyed by unique protein CDS IDs —
    // SEQKIT_SPLIT2 partitions the input so no row collides across chunks).
    // We concat with a header dedup: keep the header from the first chunk,
    // skip it on subsequent.
    //
    // Output filename matches what PHOLD_COMPARE used to emit directly, so
    // downstream consumers (MERGE_VIRAL_PHOLD, publishDir rules) don't need
    // to change.

    // Every chunk's PHOLD_COMPARE output is named phold_per_cds_predictions.tsv,
    // so `stageAs: 'per_chunk_tsv/*'` would collide on identical basenames.
    // Use `?/*` (Nextflow auto-numbered subdir, same pattern nf-core/multiqc uses)
    // to land each chunk's TSV under its own numbered subdir.
    input:
    tuple val(meta), path(per_cds_tsvs, stageAs: 'per_chunk_tsv/?/*')

    output:
    tuple val(meta), path("compare_${meta.id}/phold_per_cds_predictions.tsv"), emit: per_cds_tsv
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail
    mkdir -p compare_${meta.id}

    # Header from the first chunk, body from each. The TSV header is a single
    # line beginning with "cds_id" (per PHOLD's per_cds output schema). awk
    # NR>1 || FNR==1 keeps the first line of the first file and skips the
    # first line of every subsequent file.
    awk 'FNR > 1 || NR == 1' per_chunk_tsv/*/phold_per_cds_predictions.tsv \\
        > compare_${meta.id}/phold_per_cds_predictions.tsv

    # Row-count sanity for debugging — count is (header + sum of body rows).
    n_rows=\$(wc -l < compare_${meta.id}/phold_per_cds_predictions.tsv)
    n_chunks=\$(ls per_chunk_tsv/*/phold_per_cds_predictions.tsv | wc -l)
    echo "[merge_phold_compare] merged \${n_chunks} chunks → \${n_rows} lines (incl. header) at compare_${meta.id}/phold_per_cds_predictions.tsv" >&2

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    coreutils: \$(awk --version 2>/dev/null | head -1 | sed 's/.*GNU Awk //; s/,.*\$//')
END_VERSIONS
    """

    stub:
    """
    mkdir -p compare_${meta.id}
    printf 'cds_id\\tphrog\\tfunction\\tproduct\\tbitscore\\tannotation_confidence\\n' \\
        > compare_${meta.id}/phold_per_cds_predictions.tsv

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    coreutils: 5.3.0
END_VERSIONS
    """
}
