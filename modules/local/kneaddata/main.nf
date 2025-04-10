process KNEADDATA {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::kneaddata==0.12.0--pyhdfd78af_1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kneaddata:0.10.0--pyhdfd78af_0' :
        'quay.io/biocontainers/kneaddata:0.10.0--pyhdfd78af_0' }"

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

    input = meta.single_end ? "--input=${reads}" : "--input=${reads[0]} --input=${reads[1]}"

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

    gzip *paired_*.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kneaddata: \$(kneaddata --version | sed -e "s/kneaddata v//g")
        trimmomatic: \$(trimmomatic -version 2>&1)
        bowtie2: \$(bowtie2 --version | head -1 | cut -d' ' -f3)
    END_VERSIONS
    """
}

process PARSE_KNEADDATA {
    label 'process_small'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.22/python_3.10.sif' :
        'docker.io/biocontainers/python:3.10' }"


    input:
    tuple val(meta), path(kneaddata_log)

    output:
    tuple val(meta), path('*.csv'), emit: kneadata_stats

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    TMP_DATA=data.tmp

    grep "READ COUNT" ${kneaddata_log} | cut -d':' -f 5,7 > \$TMP_DATA

    RAW_COUNT_1="\$(grep "raw pair1" \$TMP_DATA | cut -d':' -f2 | cut -d'.' -f1 | awk ' {\$1=\$1};1')"
    #RAW_COUNT_2="\$(grep "raw pair2" \$TMP_DATA | cut -d':' -f2 | cut -d'.' -f1 | awk ' {\$1=\$1};1')"

    TRIMMED_COUNT_1="\$(grep "trimmed pair1" \$TMP_DATA | cut -d':' -f2 | cut -d'.' -f1 |  awk ' {\$1=\$1};1')"
    #TRIMMED_COUNT_2="\$(grep "trimmed pair2" \$TMP_DATA | cut -d':' -f2 | cut -d'.' -f1 | awk ' {\$1=\$1};1')"

    FINAL_COUNT_1="\$(grep "final pair1" \$TMP_DATA | cut -d':' -f2 | cut -d'.' -f1 | awk ' {\$1=\$1};1')"
    #FINAL_COUNT_2="\$(grep "final pair2" \$TMP_DATA | cut -d':' -f2 | cut -d'.' -f1 | awk ' {\$1=\$1};1')"

    echo "file,raw,trimmed,decont,final" > ${prefix}_kneaddata_stats.csv
    echo "${prefix},\$RAW_COUNT_1,\$((\$RAW_COUNT_1 - \$TRIMMED_COUNT_1)),\$((\$TRIMMED_COUNT_1 - \$FINAL_COUNT_1)),\$FINAL_COUNT_1" >> ${prefix}_kneaddata_stats.csv
    #echo "${prefix}_2,\$RAW_COUNT_2,\$((\$RAW_COUNT_2 - \$TRIMMED_COUNT_2)),\$((\$TRIMMED_COUNT_2 - \$FINAL_COUNT_2)),\$FINAL_COUNT_2" >> ${prefix}_kneaddata_stats.csv

    """
}

process SUMMARY_KNEADDATA {
    label 'process_small'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.22/python_3.10.sif' :
        'docker.io/raphsoft/python_base:3.10-R4' }"


    input:
    file(kneaddata_stats)

    output:
    path('*_mqc.html'), emit: qc_summary_plot

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    head -1 ${kneaddata_stats[0]} > kneaddata_stats.csv

    for f in \$(ls *_kneaddata_stats.csv); do
        tail -1 \$f >> kneaddata_stats.csv
    done

    qc_summary_plot.py kneaddata_stats.csv summary_report.html

    echo "<!--" > kneaddata_summary_mqc.html
    echo "id: 'kneaddata_summary'" >> kneaddata_summary_mqc.html
    echo "section_name: 'Trimming & Decontamination'" >> kneaddata_summary_mqc.html
    echo "description: 'Stats for trimmed and decontaminated reads.'" >> kneaddata_summary_mqc.html
    echo "-->" >> kneaddata_summary_mqc.html

    cat summary_report.html >> kneaddata_summary_mqc.html

    """
}


process KNEADDATA_DATABASE {
    tag "$genome_species"
    label 'process_low'

    conda "bioconda::kneaddata==0.12.0--pyhdfd78af_1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kneaddata:0.10.0--pyhdfd78af_0' :
        'quay.io/biocontainers/kneaddata:0.10.0--pyhdfd78af_0' }"

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
