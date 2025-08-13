process COVERM_MAKE {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::coverm==0.7.0--hb4818e0_2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/schirmerlab/coverm_bwamem2:0.7.0' :
        'docker.io/schirmerlab/coverm_bwamem2:0.7.0' }"

    input:
        tuple val(meta), path(reads), path(ref) // reference now can be fasta or bwa/bwamem2 index
        val ref_is_index // true if reference is a bwa/bwamem2 index

    output:
    tuple val(meta), path("*/*.bam"), emit: alignments
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def reference_stem = reference[0].getName().toString().replaceFirst(/\.[^.]+$/, '')

    input = meta.single_end ? "--single ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    out = meta.label ?: 'out'
    """
    if [[ "${ref_is_index}" == "true" ]]; then
        # reference is a bwa/bwamem2 index
        REFERENCE=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    else
        # reference is a fasta file
        REFERENCE="${ref}"
    fi

    TMPDIR=./coverm_tmp
    echo \${REFERENCE}
    coverm make $args -t $task.cpus -r \${REFERENCE} $input -o $out
    mv $out/*$prefix*.bam $out/${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CoverM: \$(coverm --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
