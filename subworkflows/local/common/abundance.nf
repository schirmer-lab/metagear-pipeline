include {BWA_INDEX} from "$projectDir/modules/nf-core/bwa/index"

include { COVERM_MAKE } from "$projectDir/modules/local/coverm/make"
include { COVERM_CONTIG_BATCH; COVERM_CONTIG_MERGE } from "$projectDir/modules/local/coverm/contig"

include { COVERM_CONTIG } from "$projectDir/modules/local/coverm/contig"



workflow ABUNDANCE {

    take:
        // label
        // reads // [meta, reads]
        // catalog // [ meta, path ]
        reads_with_index // [ meta, reads, index ]

    main:

        COVERM_MAKE ( reads_with_index, true )

        def chunkCounter = new java.util.concurrent.atomic.AtomicInteger(0)

        ch_coverm_contig = COVERM_MAKE.out.alignments
                                .map { meta, bam -> [ [id: meta.label], bam] }
                                .groupTuple (by: 0)
                                .map { meta, paths -> [ meta, paths.sort { it.toString() } ] }
                                .flatMap { meta, bams ->
                                    // Split each per-label list into fixed-size chunks (keep remainder)
                                    bams.collate(50)
                                        .withIndex()
                                        .collect { chunk, i ->
                                            def batch_id = "${meta.id}_${String.format('%03d', i+1)}"
                                            // meta carries both a unique batch id and the original label
                                            tuple( [ id: batch_id, label: meta.id, batch: i+1 ], chunk.sort { it.toString() } )  // -> [meta, [bam...]]
                                        }
                                }

        COVERM_CONTIG ( ch_coverm_contig )
        
        // helper to prepare abundance channels -> [ [id: label_suffix], [files...] ]
        def prepare_merge_channels = { suffix, ch ->
            ch
                .map { meta, file -> [ [id: meta.label], file ] }
                .groupTuple (by: 0)
                .map { tuple([id: "${it[0].id}_${suffix}"], it[1]) }
        }

        ch_coverm_merge = prepare_merge_channels( 'count', COVERM_CONTIG.out.abundance_count )
                            .concat( prepare_merge_channels( 'rpkm', COVERM_CONTIG.out.abundance_rpkm ) )
                            .concat( prepare_merge_channels( 'tpm', COVERM_CONTIG.out.abundance_tpm ) )


        COVERM_CONTIG_MERGE ( ch_coverm_merge )

        // // split merged abundance into separate channels by suffix (avoid AST/into issues)
        ch_tpm   = COVERM_CONTIG_MERGE.out.abundance_merged.filter { it[0].id.endsWith('_tpm') }
        ch_rpkm  = COVERM_CONTIG_MERGE.out.abundance_merged.filter { it[0].id.endsWith('_rpkm') }
        ch_count = COVERM_CONTIG_MERGE.out.abundance_merged.filter { it[0].id.endsWith('_count') }

        // summary channel versions
        ch_versions = COVERM_MAKE.out.versions
                        .mix(COVERM_CONTIG.out.versions)

    emit:
        alignments = COVERM_MAKE.out.alignments
        tpm = ch_tpm
        rpkm = ch_rpkm
        count = ch_count
        versions = ch_versions
}