include { CDHIT_CDHITEST } from "$projectDir/modules/local/cdhit/cdhitest/main"
include { MMSEQS_EASY_CLUSTER } from "$projectDir/modules/local/mmseqs/easy_cluster/main"

workflow CLUSTER_SEQUENCES {

    take:
        sequences
        method

    main:
        if ( method == 'cd-hit-est' ) {
            CDHIT_CDHITEST ( sequences )
            ch_clustered = CDHIT_CDHITEST.out.fasta
            ch_clusters  = CDHIT_CDHITEST.out.clusters
            ch_versions  = CDHIT_CDHITEST.out.versions

        } else if ( method == 'mmseqs2' ) {
            MMSEQS_EASY_CLUSTER ( sequences )
            ch_clustered = MMSEQS_EASY_CLUSTER.out.representatives
            ch_clusters  = MMSEQS_EASY_CLUSTER.out.clusters_tsv
            ch_versions  = MMSEQS_EASY_CLUSTER.out.versions

        } else {
            exit 1, "Unknown clustering method: ${method}"
        }

    emit:
        clustered = ch_clustered
        clusters  = ch_clusters
        versions  = ch_versions
}
