process COVERM_FILTER {
    label 'process_medium'

    conda "bioconda::coverm==0.6.1--hc216eb9_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/schirmerlab/coverm_bwamem2:0.7.0' :
        'docker.io/schirmerlab/coverm_bwamem2:0.7.0' }"

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("*/*.bam"), emit: filtered
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out = meta.label ?: 'out'
    """
    mkdir -p $out
    coverm filter $args --bam-files $bam_file --output-bam-files $out/${prefix}.bam -t $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CoverM: \$(coverm --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
