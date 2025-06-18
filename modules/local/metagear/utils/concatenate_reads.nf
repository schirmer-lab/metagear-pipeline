process CONCATENATE_READS {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path({["${meta.id}_concatenated_1.fastq.gz", "${meta.id}_concatenated_2.fastq.gz"]})

    script:
    def r1_reads = reads.findAll { it.name ==~ /.*_1\.(fq|fastq)(\.gz)?$/ }
    def r2_reads = reads.findAll { it.name ==~ /.*_2\.(fq|fastq)(\.gz)?$/ }

    def r1_out = "${meta.id}_concatenated_1.fastq.gz"
    def r2_out = "${meta.id}_concatenated_2.fastq.gz"

    // If there is more than 1 read file per id concatenate, otherwise just link to the file.
    def r1_cmd = r1_reads.size() > 1
        ? "cat ${r1_reads.join(' ')} > $r1_out"
        : "ln -s ${r1_reads[0]} $r1_out"

    def r2_cmd = r2_reads.size() > 1
        ? "cat ${r2_reads.join(' ')} > $r2_out"
        : "ln -s ${r2_reads[0]} $r2_out"

    """
    set -e
    $r1_cmd
    $r2_cmd
    """
}