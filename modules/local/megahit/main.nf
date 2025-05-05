process MEGAHIT {
    tag "$meta.id"
    label 'process_high'

    conda 'modules/nf-core/megahit/environment.yml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/megahit:1.2.9--h5b5514e_2' :
        'biocontainers/megahit:1.2.9--h5b5514e_3' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.contigs.fa.gz"), emit: contigs
    tuple val(meta), path("*.fastg")        , emit: graph
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        megahit \\
            -r ${reads} \\
            -t $task.cpus \\
            $args \\
            --out-prefix $prefix

        cat megahit_out/*.fa | sed 's/>/>${prefix}_/g' | gzip > ${prefix}.contigs.fa.gz
        rm -rf megahit_out/*.fa

        KMER_SRT="\$(ls megahit_out/intermediate_contigs/*.contigs.fa | cut -f1 -d. | rev | cut -d/ -f1 | rev | sort | uniq | head -1)"
        KMER_NUM=\${KMER_SRT#*k}

        megahit_core contig2fastg \${KMER_NUM} megahit_out/intermediate_contigs/\${KMER_SRT}.contigs.fa > ${prefix}.\${KMER_SRT}.fastg

        gzip \\
            megahit_out/intermediate_contigs/*.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    } else {
        """
        megahit \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -t $task.cpus \\
            $args \\
            --out-prefix $prefix

        cat megahit_out/*.fa | sed 's/>/>${prefix}_/g' | gzip > ${prefix}.contigs.fa.gz
        rm -rf megahit_out/*.fa

        KMER_SRT="\$(ls megahit_out/intermediate_contigs/*.contigs.fa | cut -f1 -d. | rev | cut -d/ -f1 | rev | sort | uniq | head -1)"
        KMER_NUM=\${KMER_SRT#*k}

        megahit_core contig2fastg \${KMER_NUM} megahit_out/intermediate_contigs/\${KMER_SRT}.contigs.fa > ${prefix}.\${KMER_SRT}.fastg

        gzip \\
            megahit_out/intermediate_contigs/*.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    }
}
