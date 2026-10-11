process KNEADDATA {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::kneaddata=0.12.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kneaddata:0.12.4--pyhdfd78af_0' :
        'quay.io/biocontainers/kneaddata:0.12.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path kneaddata_db

    output:
    tuple val(meta), path("*paired_{1,2}.fastq.gz"), emit: reads
    tuple val(meta), path("*kneaddata.log"), emit: kneaddata_log
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // 0.11 replaced the repeatable --input with one flag per mate.
    def input = meta.single_end ? "--unpaired ${reads}" : "--input1 ${reads[0]} --input2 ${reads[1]}"
    def db_args = kneaddata_db.collect { "-db ${it}" }.join(" ")

    """
    kneaddata \\
        ${args} \\
        --threads ${task.cpus} \\
        ${input} \\
        ${db_args} \\
        --log ${prefix}_kneaddata.log \\
        --output-prefix ${prefix} \\
        --output .

    # Named rather than globbed: KneadData keeps the contaminant reads it pulled out, and
    # *paired_* matches those too. Compressing them cost two thirds of this step on a 64 M
    # pair sample and the scratch directory then threw them away. The container has no pigz,
    # so the two that are published go in parallel instead.
    gzip -f ${prefix}_paired_1.fastq &
    gz1=\$!
    gzip -f ${prefix}_paired_2.fastq &
    gz2=\$!
    wait \$gz1
    wait \$gz2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kneaddata: \$(kneaddata --version | sed -e "s/kneaddata v//g")
        trimmomatic: \$(trimmomatic -version 2>&1)
        bowtie2: \$(bowtie2 --version | head -1 | cut -d' ' -f3)
    END_VERSIONS
    """
}

process PARSE_KNEADDATA {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.22/python_3.10.sif' :
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    tuple val(meta), path(kneaddata_log)

    output:
    tuple val(meta), path('*_read_counts.tsv'), emit: read_counts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    kneaddata_read_counts.py \\
        ${kneaddata_log} \\
        ${prefix}_read_counts.tsv \\
        --sample ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

process SUMMARY_KNEADDATA {
    label 'process_single'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.22/python_3.10.sif' :
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    path read_counts

    output:
    path "${params.workflow}_read_counts.tsv", emit: read_counts
    path "*_mqc.*"                           , emit: mqc
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    kneaddata_qc_report.py \\
        ${read_counts} \\
        --prefix ${params.workflow} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}


process KNEADDATA_DATABASE {
    tag "$genome_species"
    label 'process_low'

    conda "bioconda::kneaddata=0.12.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kneaddata:0.12.4--pyhdfd78af_0' :
        'quay.io/biocontainers/kneaddata:0.12.4--pyhdfd78af_0' }"

    input:
        tuple val(genome_species), val(dbtype)

    output:
        tuple val(genome_species), path("$genome_species"), emit: database
        path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    kneaddata_database --download $genome_species $dbtype "$genome_species"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kneaddata: \$(kneaddata --version | sed -e "s/kneaddata v//g")
    END_VERSIONS
    """
}
