process MMSEQS_EASY_CLUSTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:17.b804f--hd6d6fdc_1':
        'biocontainers/mmseqs2:17.b804f--hd6d6fdc_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.clusters.tsv"), emit: clusters_tsv
    tuple val(meta), path("*.representative.fa.gz"), emit: representatives
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}

    if [[ $fasta == *.gz ]]
    then
        OUT_NAME=\$(basename ${fasta})
        OUT_NAME=\${OUT_NAME%.gz}
    fi
    OUT_NAME=\${OUT_NAME%.*}


    mmseqs easy-cluster ${fasta} cluster_res ${prefix} \\
        $args \\
        --threads ${task.cpus}

    mv cluster_res_cluster.tsv \${OUT_NAME}.clusters.tsv
    mv cluster_res_rep_seq.fasta \${OUT_NAME}.representative.fa
    mv cluster_res_all_seqs.fasta \${OUT_NAME}.all_seqs.fa

    gzip \${OUT_NAME}.representative.fa \${OUT_NAME}.all_seqs.fa


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    touch ${prefix}/${prefix}.{0..9}
    touch ${prefix}/${prefix}.dbtype
    touch ${prefix}/${prefix}.index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
