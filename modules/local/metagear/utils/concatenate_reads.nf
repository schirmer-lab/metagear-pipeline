process CONCATENATE_READS {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path({"${meta.id}_concatenated_[12].fastq.gz"})

    script:
    def r1_reads = reads.findAll { it.name ==~ /.*_(1|R1_\d+)\.(fq|fastq)(\.gz)?$/ }
    def r2_reads = reads.findAll { it.name ==~ /.*_(2|R2_\d+)\.(fq|fastq)(\.gz)?$/ }

    def r1_out = "${meta.id}_concatenated_1.fastq.gz"
    def r1_cmd = r1_reads.size() > 1
        ? "cat ${r1_reads.join(' ')} > $r1_out"
        : "ln -s ${r1_reads[0]} $r1_out"

    // Prepare r2 only if available
    def r2_cmd = ""
    def r2_out = ""
    if (r2_reads) {
        r2_out = "${meta.id}_concatenated_2.fastq.gz"
        r2_cmd = r2_reads.size() > 1
            ? "cat ${r2_reads.join(' ')} > $r2_out"
            : "ln -s ${r2_reads[0]} $r2_out"
    }

    """
    set -e
    $r1_cmd
    $r2_cmd
    """
}
