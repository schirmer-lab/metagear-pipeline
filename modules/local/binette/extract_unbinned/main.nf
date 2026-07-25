process EXTRACT_UNBINNED {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::seqkit=2.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.2.0--h9ee0642_0':
        'biocontainers/seqkit:2.2.0--h9ee0642_0' }"

    input:
    // chromosome_fa  — the per-sample chromosome partition fed into binning
    // bin_fastas     — list of Binette's final_bins/*.fa (the MIMAG-MQ+ set)
    // When bin_fastas is empty (sample produced zero passing bins), every
    // chromosome contig is unbinned and we emit the chromosome FASTA as-is.
    tuple val(meta), path(chromosome_fa), path(bin_fastas, stageAs: 'bin_fastas/*')

    output:
    tuple val(meta), path("${prefix}.unbinned.fna.gz"), emit: unbinned
    tuple val(meta), path("${prefix}.binned_ids.txt") , emit: binned_ids
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    # Collect every contig header that ended up in a Binette final bin.
    # Bin FASTAs are plain .fa (Binette default); handle either plain or gzipped just in case.
    # Note: Nextflow runs .command.sh with `bash -C` (noclobber) so we must NOT
    # pre-create the output file with `: > <file>` — the pipeline below redirects
    # into it directly, and a pre-existing file would trigger
    # "cannot overwrite existing file".
    shopt -s nullglob
    for f in bin_fastas/*.fa bin_fastas/*.fa.gz bin_fastas/*.fasta bin_fastas/*.fasta.gz; do
        if [[ "\$f" == *.gz ]]; then
            zcat "\$f"
        else
            cat "\$f"
        fi
    done | awk '/^>/ {sub(/^>/, ""); split(\$0, a, /[ \\t]/); print a[1]}' \\
         | sort -u >| ${prefix}.binned_ids.txt

    if [ -s ${prefix}.binned_ids.txt ]; then
        seqkit grep -v -n -f ${prefix}.binned_ids.txt $args ${chromosome_fa} \\
            | gzip --no-name > ${prefix}.unbinned.fna.gz
    else
        # Sample produced zero bins — every chromosome contig is unbinned.
        cp ${chromosome_fa} ${prefix}.unbinned.fna.gz
    fi

    # If the unbinned FASTA is empty (every contig binned), drop the file so
    # downstream `optional: true` consumers see no item rather than an empty FASTA.
    if [ ! -s ${prefix}.unbinned.fna.gz ] || [ "\$(zcat ${prefix}.unbinned.fna.gz | head -c1 | wc -c)" = "0" ]; then
        rm -f ${prefix}.unbinned.fna.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.binned_ids.txt
    echo "" | gzip > ${prefix}.unbinned.fna.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: 2.2.0
    END_VERSIONS
    """
}
