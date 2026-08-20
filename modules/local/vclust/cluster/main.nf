process VCLUST_CLUSTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vclust:1.3.1--py313h9ee0642_0':
        'biocontainers/vclust:1.3.1--py313h9ee0642_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.clusters.tsv"), emit: clusters_tsv
    tuple val(meta), path("*.representative.fa.gz"), emit: representatives
    tuple val(meta), path("*.ani.tsv"), emit: ani
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefilter_args = task.ext.args   ?: '--min-ident 0.90'
    def align_args     = task.ext.args2  ?: '--outfmt lite'
    def cluster_args   = task.ext.args3  ?: '--algorithm single --metric ani --ani 0.95 --rcov 0.85'
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    INPUT=${fasta}
    OUT_NAME=\$(basename ${fasta})
    if [[ ${fasta} == *.gz ]]
    then
        OUT_NAME=\${OUT_NAME%.gz}
        gunzip -c ${fasta} > ${prefix}.fna
        INPUT=${prefix}.fna
    fi
    OUT_NAME=\${OUT_NAME%.*}

    # Three stages, as Vclust requires. The prefilter screens candidate pairs by
    # k-mer identity at a threshold deliberately below the clustering threshold,
    # so no pair that could pass clustering is discarded before alignment.
    vclust prefilter -i \${INPUT} -o \${OUT_NAME}.prefilter.txt \\
        ${prefilter_args} \\
        -t ${task.cpus}

    vclust align -i \${INPUT} -o \${OUT_NAME}.ani.tsv \\
        --filter \${OUT_NAME}.prefilter.txt \\
        ${align_args} \\
        -t ${task.cpus}

    vclust cluster -i \${OUT_NAME}.ani.tsv -o \${OUT_NAME}.vclust.tsv \\
        --ids \${OUT_NAME}.ani.ids.tsv \\
        ${cluster_args}

    # Vclust emits <sequence>\\t<cluster_id> keyed by a numeric cluster. Every
    # consumer downstream of CLUSTER_SEQUENCES expects MMseqs2 easy-cluster's
    # shape instead (representative<TAB>member, plus a representative FASTA), so
    # the conversion happens here and nothing downstream needs to change.
    vclust_to_pairs.py \\
        --fasta \${INPUT} \\
        --clusters \${OUT_NAME}.vclust.tsv \\
        --out-pairs \${OUT_NAME}.clusters.tsv \\
        --out-representatives \${OUT_NAME}.representative.fa.gz

    if [[ ${fasta} == *.gz ]]
    then
        rm ${prefix}.fna
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vclust: \$(vclust --version 2>&1 | sed 's/^v//')
    END_VERSIONS
    """
}
