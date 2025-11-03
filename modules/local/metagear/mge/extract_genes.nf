process EXTRACT_GENES {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqtk=1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/schirmerlab/python310:25.09.10' :
        'docker.io/schirmerlab/python310:25.09.10' }"

    input:
    tuple val(meta), path(contig_ids), path(all_genes)

    output:
    tuple val(meta), path("*.fasta"), emit: extracted_genes, optional: true
    tuple val(meta), path("*.ids.txt"), emit: extracted_gene_ids, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    PREFIX=${prefix} extract_genes.sh ${contig_ids} ${all_genes}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk 2>&1 | awk '/Version:/ {print \$2; exit}')
    END_VERSIONS
    """
}

