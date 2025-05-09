process METAPHLAN_METAPHLAN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metaphlan:4.1.1--pyhdfd78af_0' :
        'biocontainers/metaphlan:4.1.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input)
    path metaphlan_db_latest

    output:
    tuple val(meta), path("*_microbial_profile.txt"), emit: microbial_profile
    tuple val(meta), path("*_viral_profile.txt"), optional: true, emit:viral_profile
    tuple val(meta), path("*.biom"), emit: biom
    tuple val(meta), path('*.bowtie2out.txt'), optional:true, emit: bt2out
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_type = "$input" =~ /.*\.(fastq|fq)/ ? "--input_type fastq" : "$input" =~ /.*\.(fasta|fna|fa)/ ? "--input_type fasta" : "$input".endsWith(".bowtie2out.txt") ? "--input_type bowtie2out" : "--input_type sam"
    def input_data  = ("$input_type".contains("fastq")) && !meta.single_end ? "${input[0]},${input[1]}" : "$input"
    def bowtie2_out = "$input_type" == "--input_type bowtie2out" || "$input_type" == "--input_type sam" ? '' : "--bowtie2out ${prefix}.bowtie2out.txt"

    // def profile_virus =  --profile_vsc
    def virus  = "$args".contains("--profile_vsc") ? "--profile_vsc --vsc_out ${prefix}_viral_profile.txt" : ''

    """
    BT2_DB=`find -L "${metaphlan_db_latest}" -name "*rev.1.bt2*" -exec dirname {} \\;`
    BT2_DB_INDEX=`find -L ${metaphlan_db_latest} -name "*.rev.1.bt2*" | sed 's/\\.rev.1.bt2.*\$//' | sed 's/.*\\///'`

    metaphlan \\
        --nproc $task.cpus \\
        $input_type \\
        $input_data \\
        $args \\
        $virus \\
        $bowtie2_out \\
        --bowtie2db \$BT2_DB \\
        --index \$BT2_DB_INDEX \\
        --biom ${prefix}.biom \\
        --output_file ${prefix}_microbial_profile.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: \$(metaphlan --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """
}

process METAPHLAN_MERGE_PROFILES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metaphlan:4.1.1--pyhdfd78af_0' :
        'biocontainers/metaphlan:4.1.1--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(profiles)

    output:
        tuple val(meta), path("merged_*_profiles.txt") , emit: merged_profiles
        path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    if [ "$prefix" = "microbial" ]; then
        merge_metaphlan_tables.py $profiles > merged_microbial_profiles.txt
    fi

    if [ "$prefix" = "viral" ]; then
        merge_vsc_tables.py $profiles > merged_viral_profiles.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: \$(metaphlan --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """
}
