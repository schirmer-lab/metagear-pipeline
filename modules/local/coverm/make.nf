process COVERM_MAKE {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::coverm==0.7.0--hb4818e0_2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/schirmerlab/coverm_bwamem2:0.7.0' :
        'docker.io/schirmerlab/coverm_bwamem2:0.7.0' }"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/coverm:0.7.0--hb4818e0_2' :
    //     'biocontainers/coverm:0.7.0--hb4818e0_2' }"

    input:
        tuple val(meta), path(reads), path(reference) // Genome now can be fasta or bwa/bwamem2 index


    output:
    tuple val(meta), path("*/*.bam"), emit: alignments
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference_stem = reference[0].getName().toString().replaceFirst(/\.[^.]+$/, '')

    input = meta.single_end ? "--single ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    out = meta.label ?: 'out'
    """
    TMPDIR=./coverm_tmp
    echo ${reference_stem}
    coverm make $args -t $task.cpus -r ${reference_stem} $input -o $out
    mv $out/*$prefix*.bam $out/${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CoverM: \$(coverm --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
