process COVERM_CONTIG {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::coverm==0.7.0--hb4818e0_2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coverm:0.7.0--hb4818e0_2' :
        'biocontainers/coverm:0.7.0--hb4818e0_2' }"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.abundance_count.tsv"), emit: abundance_count
    tuple val(meta), path("*.abundance_trimmed_mean.tsv"), emit: abundance_trimmed_mean
    tuple val(meta), path("*.abundance_rpkm.tsv"), emit: abundance_rpkm
    tuple val(meta), path("*.abundance_tpm.tsv"), emit: abundance_tpm
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    coverm contig --methods count --bam-files $bams -t $task.cpus $args 1> ${prefix}.abundance_count.tsv 2> log_count.txt
    sed -i '1 s/ Read Count//g' ${prefix}.abundance_count.tsv

    coverm contig --methods trimmed_mean --bam-files $bams -t $task.cpus $args2 1> ${prefix}.abundance_trimmed_mean.tsv 2> log_trimmed_mean.txt
    sed -i '1 s/ Trimmed Mean//g' ${prefix}.abundance_trimmed_mean.tsv

    #coverm contig --methods trimmed_mean --bam-files $bams -t $task.cpus 1> ${prefix}.abundance_trimmed_mean_unfiltered.tsv 2> log_trimmed_mean_unfiltered.txt
    #sed -i '1 s/ Trimmed Mean//g' ${prefix}.abundance_trimmed_mean_unfiltered.tsv

    coverm contig --methods rpkm --bam-files $bams -t $task.cpus $args2 1> ${prefix}.abundance_rpkm.tsv 2> log_rpkm.txt
    sed -i '1 s/ RPKM//g' ${prefix}.abundance_rpkm.tsv

    #coverm contig --methods rpkm --bam-files $bams -t $task.cpus 1> ${prefix}.abundance_rpkm_unfiltered.tsv 2> log_rpkm_unfiltered.txt
    #sed -i '1 s/ RPKM//g' ${prefix}.abundance_rpkm_unfiltered.tsv

    coverm contig --methods tpm --bam-files $bams -t $task.cpus $args2 1> ${prefix}.abundance_tpm.tsv 2> log_tpm.txt
    sed -i '1 s/ TPM//g' ${prefix}.abundance_tpm.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CoverM: \$(coverm --version | cut -d' ' -f2)
    END_VERSIONS
    """
}

process COVERM_CONTIG_BATCH {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::coverm==0.7.0--hb4818e0_2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coverm:0.7.0--hb4818e0_2' :
        'biocontainers/coverm:0.7.0--hb4818e0_2' }"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.abundance_count.tsv"), emit: abundance_count
    tuple val(meta), path("*.abundance_rpkm.tsv"), emit: abundance_rpkm
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """

    files=($bams)

    mkdir tmp_batches

    for ((i=0; i<\${#files[@]}; i+=50)); do
        batch_files="\${files[@]:i:50}"
        batch_prefix=${prefix}_batch_\${i}

        coverm contig --methods count --bam-files \$batch_files -t $task.cpus $args 1> \${batch_prefix}.abundance_count.tsv 2> \${batch_prefix}_log_count.txt
        sed -i '1 s/ Read Count//g' \${batch_prefix}.abundance_count.tsv

        coverm contig --methods rpkm --bam-files \$batch_files -t $task.cpus $args2 1> \${batch_prefix}.abundance_rpkm.tsv 2> \${batch_prefix}_log_rpkm.txt
        sed -i '1 s/ RPKM//g' \${batch_prefix}.abundance_rpkm.tsv

    done

    # Remove the temporary files
    rm -rf tmp_batches

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CoverM: \$(coverm --version | cut -d' ' -f2)
    END_VERSIONS
    """
}


process COVERM_CONTIG_MERGE {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.22/python_3.10.sif' :
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    tuple val(meta), path(tsv_files)

    output:
    tuple val(meta), path("*_merged.tsv"), emit: abundance_merged
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    coverm_merge.py ${tsv_files} -o ${prefix}_merged.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
