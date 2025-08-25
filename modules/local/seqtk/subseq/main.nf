process SEQTK_SUBSEQ {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::seqtk=1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1' :
        'biocontainers/seqtk:1.4--he4a0461_1' }"

    input:
    tuple val(meta), path(sequences), path(filter_list)

    output:
    tuple val(meta), path('*.gz'), optional: true, emit: sequences
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def args2 = task.ext.args2 ?: '.default' // Name for the resulting contig
    def m = args2 =~ /--catalog_name\s+(\S+)/
    def catalog_name = m.find() ? m.group(1) : ''

    def ext = "fa"
    if ("$sequences" ==~ /.+\.fq|.+\.fq.gz|.+\.fastq|.+\.fastq.gz/) {
        ext = "fq"
    }
    """
    zcat $sequences | \\
    seqtk \\
        subseq \\
        $args \\
        - \\
        $filter_list > sequences.fa

    if [ -s sequences.fa ]; then
        cat sequences.fa | gzip --no-name > "${prefix}.${catalog_name}.gz"
        rm -rf sequences.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}

process SEQTK_SPLIT_BY_LENGTH {
    tag "$sequences"
    label 'process_single'

    conda "bioconda::seqkit=2.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.2.0--h9ee0642_0':
        'biocontainers/seqkit:2.2.0--h9ee0642_0' }"
    input:
    tuple val(meta), path(sequences)
    val min_contig_length

    output:
    tuple val(meta), path('*_long.fasta.gz'), optional: true, emit: long_sequences
    tuple val(meta), path('*_short.fasta.gz'), optional: true, emit: short_sequences
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    INPUT=$sequences
    if [[ $sequences == *.gz ]]
    then
        INPUT=\$(basename ${sequences})
        INPUT=\${INPUT%.gz}
        gunzip -c $sequences > \$INPUT
    fi

    # Long contigs (>= min_contig_length bp)
    seqkit seq $args -m ${min_contig_length}  \${INPUT} | gzip > \${INPUT%.*}_long.fasta.gz

    # Short contigs (< min_contig_length bp)
    seqkit seq $args -M ${min_contig_length} \${INPUT} | gzip > \${INPUT%.*}_short.fasta.gz

    rm -rf \$INPUT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}

