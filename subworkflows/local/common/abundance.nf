include {BWA_INDEX} from "$projectDir/modules/nf-core/bwa/index"

include { COVERM_MAKE } from "$projectDir/modules/local/coverm/make"
include { COVERM_CONTIG_BATCH; COVERM_CONTIG_BATCH_MERGE as COVERM_BATCH_MERGE_COUNT; COVERM_CONTIG_BATCH_MERGE as COVERM_BATCH_MERGE_RPKM } from "$projectDir/modules/local/coverm/contig"


workflow ABUNDANCE {

    take:
        label
        reads // [meta, reads]
        catalog // [ meta, path ]

    main:
        BWA_INDEX ( catalog ) // build index

        ch_index = BWA_INDEX.out.index.map {[it[1]]}
        ch_make = reads.combine(ch_index)

        COVERM_MAKE ( ch_make, true )

        COVERM_MAKE.out.alignments
            .map { it -> it[1] }
            .collect()
            .map { it -> tuple([id: label], it)}
            .set { ch_coverm_contig }

        // generate summary table, e.g. count, rpkm, tpm
        COVERM_CONTIG_BATCH ( ch_coverm_contig )

        ch_coverm_contig_merge_count = COVERM_CONTIG_BATCH.out.abundance_count.map { it -> tuple( [id: it[0].id + '_count'], it[1] ) }
        COVERM_BATCH_MERGE_COUNT ( ch_coverm_contig_merge_count )

        ch_coverm_contig_merge_rpkm = COVERM_CONTIG_BATCH.out.abundance_rpkm.map { it -> tuple( [id: it[0].id + '_rpkm'], it[1] ) }
        COVERM_BATCH_MERGE_RPKM ( ch_coverm_contig_merge_rpkm)

        // summary channel versions
        ch_versions = BWA_INDEX.out.versions
                        .mix(COVERM_MAKE.out.versions)
                        .mix(COVERM_CONTIG_BATCH.out.versions)

    emit:
        catalog_index = BWA_INDEX.out.index
        alignments = COVERM_MAKE.out.alignments
        abundance_count = COVERM_BATCH_MERGE_COUNT.out.abundance_count
        abundance_rpkm = COVERM_BATCH_MERGE_RPKM.out.abundance_count
        versions = ch_versions
}

