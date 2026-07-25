process MERGE_CONTIG_CLASSIFICATION {
    tag "$meta.id"
    label 'process_single'

    // Pure-python helper (bin/merge_contig_classification.py); container only
    // needs Python 3 + the stdlib.
    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11':
        'biocontainers/python:3.11' }"

    input:
    // All evidence channels joined upstream into one tuple per sample:
    //   contigs       — chromosome FASTA (defines the universe of contig_ids)
    //   viral_ids     — MERGE_TABLES viral_ids_to_keep.txt (may be empty file)
    //   plasmid_ids   — MERGE_TABLES plasmid_ids_to_keep.txt (may be empty file)
    //   bin_fastas    — Binette MAG FASTAs (list of paths) staged as bins/<name>.fa
    //                   so the python script can iterate `bins/` like a real
    //                   directory. When the sample has no bins, the upstream
    //                   join falls back to assets/empty.txt (a single
    //                   non-FASTA file); it stages as bins/empty.txt and the
    //                   python script's .fa-extension filter skips it.
    //   mmseqs_lca    — MMSEQS_EASYTAXONOMY ${prefix}_lca.tsv (may be missing if no unbinned)
    //   tiara         — TIARA_TIARA ${prefix}.txt (may be missing if --chromosome_dir bypass)
    //
    // Each "missing" channel is staged as the assets/empty.txt fallback (see
    // the subworkflow's `.join(remainder:true)` plumbing); the python
    // script silently treats empty/missing-dir bin sets + zero-byte ID
    // files as empty evidence sets.
    tuple val(meta), path(contigs), path(viral_ids), path(plasmid_ids), path(bin_fastas, stageAs: 'bins/*'), path(mmseqs_lca), path(tiara)

    output:
    tuple val(meta), path("${prefix}.contigs.tsv"), emit: tsv
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    // Resolve optional inputs to either the staged path or omit the flag.
    // assets/empty.txt is the wrapper for "this evidence channel was empty";
    // the script treats zero-byte input files as empty sets, so we can pass
    // them through unconditionally.
    """
    merge_contig_classification.py \\
        --contigs ${contigs} \\
        --sample "${meta.id}" \\
        --viral-ids ${viral_ids} \\
        --plasmid-ids ${plasmid_ids} \\
        --bin-dir bins \\
        --mmseqs-lca ${mmseqs_lca} \\
        --tiara ${tiara} \\
        --output ${prefix}.contigs.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'contig_id\\tsample\\tlength\\tprimary_class\\tclassifier\\tlineage\\tconfidence\\tbin_id\\tgenomad_viral\\tgenomad_plasmid\\ttiara_label\\n' > ${prefix}.contigs.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
