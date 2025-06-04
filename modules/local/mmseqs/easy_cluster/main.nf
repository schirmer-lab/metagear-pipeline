process MMSEQS_EASY_CLUSTER {
   
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.clusters.tsv"), emit: clusters_tsv
    tuple val(meta), path("*.representative.fa.gz"), emit: representatives
    path "versions.yml" , emit: versions

    """
    echo "test" > test.clusters.tsv
    echo "test" > test.representative.fa.gz
    echo "test" > versions.yml
    """
}
