process BINETTE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::binette=1.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/binette:1.2.1--pyh106432d_1':
        'quay.io/biocontainers/binette:1.2.1--pyh106432d_1' }"

    input:
    // Binette 1.2.1's CLI requires --bin_dirs XOR --contig2bin_tables (not both).
    // Both upstream binners (SemiBin2 + MetaBAT2) emit gzipped bin FASTAs in
    // directories, so we stage each into its own subdir and pass them as a
    // space-separated list to --bin_dirs.
    tuple val(meta), path(contigs), path(semibin_bins, stageAs: 'semibin_bins/*'), path(metabat_bins, stageAs: 'metabat_bins/*')
    path checkm2_db

    output:
    tuple val(meta), path("${prefix}/final_bins/*.fa")                       , emit: bins         , optional: true
    tuple val(meta), path("${prefix}/final_bins_quality_reports.tsv")        , emit: quality
    tuple val(meta), path("${prefix}/final_contig_to_bin.tsv")               , emit: contig_to_bin
    path "versions.yml"                                                       , emit: versions
    // Note: both directory emits (`outdir` for the whole ${prefix}/ and
    // `bins_dir` for ${prefix}/final_bins/) were dropped. They shadowed the
    // file-glob bins emit in publishDir — Nextflow processes the directory
    // emit's view of the same files and silently skips them when the dir's
    // saveAs returns null. With only file-typed emits (matches the
    // FIND_REPRESENTATIVES pattern), publishDir resolves cleanly per file.
    //
    // Downstream consumers that previously took bins_dir as a directory now
    // receive the bins file list and stage it as a directory via
    // `path(bins, stageAs: 'bins/*')` — the consumer's view is unchanged.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    # Build --bin_dirs from whichever staged directories actually contain bin
    # FASTAs. When SemiBin or MetaBAT2 produced 0 bins for this sample,
    # bacterial_binning.nf passes an empty list (via remainder:true + ?:[]);
    # the staged directory is either absent or empty, and nullglob makes the
    # *.fa* expansion safely empty. Same pattern as EXTRACT_UNBINNED's
    # bin_fastas/ handling.
    shopt -s nullglob
    semibin_files=( semibin_bins/*.fa semibin_bins/*.fa.gz semibin_bins/*.fasta semibin_bins/*.fasta.gz )
    metabat_files=( metabat_bins/*.fa metabat_bins/*.fa.gz metabat_bins/*.fasta metabat_bins/*.fasta.gz )

    bin_dirs_arg=""
    [ \${#semibin_files[@]} -gt 0 ] && bin_dirs_arg="\$bin_dirs_arg semibin_bins"
    [ \${#metabat_files[@]} -gt 0 ] && bin_dirs_arg="\$bin_dirs_arg metabat_bins"

    # Edge case: both binners produced 0 bins for this sample (very low
    # complexity input). Emit empty deliverables and skip Binette — its
    # `bins` output is already declared optional, and classification.nf:174
    # handles "no bins" downstream via remainder:true (EXTRACT_UNBINNED
    # treats every chromosome contig as unbinned in that case).
    if [ -z "\${bin_dirs_arg// /}" ]; then
        echo "[BINETTE] Both binners produced 0 bins for ${meta.id} — emitting empty reports" >&2
        mkdir -p ${prefix}/final_bins
        : > ${prefix}/final_bins_quality_reports.tsv
        : > ${prefix}/final_contig_to_bin.tsv
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            binette: skipped_no_input_bins
        END_VERSIONS
        exit 0
    fi

    binette \\
        --contigs ${contigs} \\
        --bin_dirs \$bin_dirs_arg \\
        --checkm2_db ${checkm2_db} \\
        --outdir ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: \$(binette --version 2>&1 | sed 's/.*[Bb]inette //; s/[,\\s].*//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/final_bins
    touch ${prefix}/final_bins/bin_1.fa
    touch ${prefix}/final_bins_quality_reports.tsv
    touch ${prefix}/final_contig_to_bin.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: 1.2.1
    END_VERSIONS
    """
}
