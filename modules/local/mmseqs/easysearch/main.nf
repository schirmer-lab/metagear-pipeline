process MMSEQS_EASYSEARCH {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mmseqs2:17.b804f--hd6d6fdc_1'
        : 'biocontainers/mmseqs2:17.b804f--hd6d6fdc_1'}"

    input:
    tuple val(meta), path(fasta), path(db_target)

    output:
    tuple val(meta), path("${prefix}_search.tsv"), emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // def args2 = task.ext.args2 ?: "*.dbtype"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    # Extract files with specified args based suffix | remove suffix | isolate longest common substring of files

    INPUT=$db_target
    if [[ $db_target == *.gz ]]
    then
        gunzip -c $db_target >| \$PWD/${prefix}_target.fa
        DB_TARGET=\$PWD/${prefix}_target.fa
    fi

    mmseqs \\
        easy-search \\
        ${fasta} \\
        \${DB_TARGET} \\
        ${prefix}_search.tsv \\
        tmp1 \\
        ${args} \\
        --threads ${task.cpus}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
