process GENOMAD_EXTRACT_CHROMOSOME {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::seqkit=2.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.2.0--h9ee0642_0':
        'biocontainers/seqkit:2.2.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(contigs), path(viral_ids), path(plasmid_ids)

    output:
    tuple val(meta), path("*.chromosome.fna.gz"), emit: chromosome_fasta, optional: true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${viral_ids} ${plasmid_ids} | sort -u > exclude_ids.txt

    if [ -s exclude_ids.txt ]; then
        seqkit grep -v -n -f exclude_ids.txt $args ${contigs} \\
            | gzip --no-name > "${prefix}.chromosome.fna.gz"
    else
        cp ${contigs} "${prefix}.chromosome.fna.gz"
    fi

    if [ ! -s "${prefix}.chromosome.fna.gz" ]; then
        rm -f "${prefix}.chromosome.fna.gz"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > "${prefix}.chromosome.fna.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: 2.2.0
    END_VERSIONS
    """
}
