process MMSEQS_TAXONOMY {
    tag "$meta.id"
    container "/nfs/data/database/singularity/denglab-viroprofiler-base-v0.2.img"

    input:
    tuple val(meta), path(contigs)
    path mmseqs_db

    output:
    path "*.mmseqsTaxaRst.tsv", emit: taxa_mmseqs_ch
    path "*.mmseqsTaxaRst_report.*"
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mmseqs createdb $contigs qry
    mmseqs taxonomy qry ${mmseqs_db}/refseq_viral mmseqsTaxaRst tmp --tax-lineage 1 --majority 0.4 --vote-mode 1 --lca-mode 3 --orf-filter 0 --threads $task.cpus

    mmseqs createtsv qry mmseqsTaxaRst ${prefix}.mmseqsTaxaRst.tsv
    mmseqs taxonomyreport ${mmseqs_db}/refseq_viral mmseqsTaxaRst ${prefix}.mmseqsTaxaRst_report.txt --report-mode 0
    mmseqs taxonomyreport ${mmseqs_db}/refseq_viral mmseqsTaxaRst ${prefix}.mmseqsTaxaRst_report.html --report-mode 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MMseqs2: \$(grep "MMseqs Version" .command.log | head -n1 | sed 's/.*\t//g')
    END_VERSIONS
    """
}
