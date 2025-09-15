include { VAMB_CONCATENATE_FASTA } from "$projectDir/modules/local/vamb/main"

include { CDHIT_CDHITEST } from "$projectDir/modules/local/cdhit/cdhitest/main"
include { MMSEQS_EASY_CLUSTER } from "$projectDir/modules/local/mmseqs/easy_cluster/main"

workflow CLUSTER_SEQUENCES {

    take:
        sequences // tuple (meta, fasta) -> [ [id: ab12, label: ABD ], fasta ] :: label is REQUIRED for clustering
        method
        concatenate

    main:

        ch_sequences = sequences

        if ( concatenate ) {

            ch_merged_genes = sequences.map{ [[ id: it[0].label ?: it[0].id ], it[1] ] }
                .groupTuple(by: 0)
                .map { meta, paths -> [ meta, paths.sort { it.toString() } ] }

            VAMB_CONCATENATE_FASTA ( ch_merged_genes )
            ch_sequences = VAMB_CONCATENATE_FASTA.out.catalog

        }

        if ( method == 'cd-hit-est' ) {
            CDHIT_CDHITEST ( ch_sequences )
            ch_clustered = CDHIT_CDHITEST.out.fasta
            ch_clusters  = CDHIT_CDHITEST.out.clusters
            ch_versions  = CDHIT_CDHITEST.out.versions

        } else if ( method == 'mmseqs2' ) {
            MMSEQS_EASY_CLUSTER ( ch_sequences )
            ch_clustered = MMSEQS_EASY_CLUSTER.out.representatives
            ch_clusters  = MMSEQS_EASY_CLUSTER.out.clusters_tsv
            ch_versions  = MMSEQS_EASY_CLUSTER.out.versions

        } else {
            exit 1, "Unknown clustering method: ${method}"
        }

    emit:
        representative = ch_clustered
        clusters  = ch_clusters
        versions  = ch_versions
}
