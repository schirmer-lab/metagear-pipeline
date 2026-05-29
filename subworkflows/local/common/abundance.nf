include { BWA_INDEX } from "$projectDir/modules/nf-core/bwa/index"

include { COVERM_MAKE } from "$projectDir/modules/local/coverm/make"
include { COVERM_CONTIG_BATCH; COVERM_CONTIG_MERGE } from "$projectDir/modules/local/coverm/contig"

include { COVERM_CONTIG } from "$projectDir/modules/local/coverm/contig"
include { COVERM_GENOME } from "$projectDir/modules/local/coverm/genome"


workflow ABUNDANCE {

    take:
        reads_with_sequences // [ meta, reads, catalog ]: meta must contain id and label
        mode                 // 'contig' | 'genome' — selects coverm subcommand
        genome_definition    // path: contig→genome TSV consumed when mode == 'genome'.
                             // In contig mode pass a placeholder (e.g.
                             // file("$projectDir/assets/empty.txt")); it is never read.

    main:

        /* -- Build index for abundance estimation (only once) --- */
        ch_catalogs = reads_with_sequences
                        .map { meta, reads, catalog -> [ [id: meta.label, src: catalog.getName()], catalog ] }
                        .unique()

        BWA_INDEX ( ch_catalogs )

        // Combine back the index with reads and catalog for abundance estimation
        ch_reads_with_index = reads_with_sequences
                        .map { meta, reads, catalog -> [ [id: meta.label, src: catalog.getName()], meta, reads ] }
                        .combine( BWA_INDEX.out.index, by: 0 )
                        .map { src, meta, reads, index -> [ meta, reads, index ] }

        COVERM_MAKE ( ch_reads_with_index, true )

        def chunkCounter = new java.util.concurrent.atomic.AtomicInteger(0)

        ch_coverm_contig = COVERM_MAKE.out.alignments
                                .map { meta, bam -> [ [id: meta.label], bam] }
                                .groupTuple (by: 0)
                                .map { meta, paths -> [ meta, paths.sort { it.toString() } ] }
                                .flatMap { meta, bams ->
                                    // Split each per-label list into fixed-size chunks (keep remainder)
                                    bams.collate(params.files_batch_size)
                                        .withIndex()
                                        .collect { chunk, i ->
                                            def batch_id = "${meta.id}_${String.format('%03d', i+1)}"
                                            // meta carries both a unique batch id and the original label
                                            tuple( [ id: batch_id, label: meta.id, batch: i+1 ], chunk.sort { it.toString() } )  // -> [meta, [bam...]]
                                        }
                                }

        // Branch on mode. The contig arm is byte-identical to the pre-mode flow
        // (same process, same input channel) so existing genes /
        // virus / msp -resume runs cache-hit.
        if ( mode == 'genome' ) {
            // Attach the cohort-global contig→genome definition to every batch.
            // `combine` against a value/queue channel of one item broadcasts.
            COVERM_GENOME ( ch_coverm_contig.combine( genome_definition ) )
            ch_abund_count = COVERM_GENOME.out.abundance_count
            ch_abund_rpkm  = COVERM_GENOME.out.abundance_rpkm
            ch_abund_tpm   = COVERM_GENOME.out.abundance_tpm
            ch_abund_versions = COVERM_GENOME.out.versions
        } else {
            COVERM_CONTIG ( ch_coverm_contig )
            ch_abund_count = COVERM_CONTIG.out.abundance_count
            ch_abund_rpkm  = COVERM_CONTIG.out.abundance_rpkm
            ch_abund_tpm   = COVERM_CONTIG.out.abundance_tpm
            ch_abund_versions = COVERM_CONTIG.out.versions
        }

        // helper to prepare abundance channels -> [ [id: label_suffix], [files...] ]
        def prepare_merge_channels = { suffix, ch ->
            ch
                .map { meta, file -> [ [id: meta.label], file ] }
                .groupTuple (by: 0)
                .map { tuple([id: "${it[0].id}_${suffix}"], it[1]) }
        }

        ch_coverm_merge = prepare_merge_channels( 'count', ch_abund_count )
                            .concat( prepare_merge_channels( 'rpkm', ch_abund_rpkm ) )
                            .concat( prepare_merge_channels( 'tpm', ch_abund_tpm ) )


        COVERM_CONTIG_MERGE ( ch_coverm_merge )

        // // split merged abundance into separate channels by suffix (avoid AST/into issues)
        ch_tpm   = COVERM_CONTIG_MERGE.out.abundance_merged.filter { it[0].id.endsWith('_tpm') }
        ch_rpkm  = COVERM_CONTIG_MERGE.out.abundance_merged.filter { it[0].id.endsWith('_rpkm') }
        ch_count = COVERM_CONTIG_MERGE.out.abundance_merged.filter { it[0].id.endsWith('_count') }

        // summary channel versions
        ch_versions = COVERM_MAKE.out.versions
                        .mix(ch_abund_versions)
        // ch_versions = Channel.empty()

    emit:
        index = BWA_INDEX.out.index
        alignments = COVERM_MAKE.out.alignments
        tpm = ch_tpm
        rpkm = ch_rpkm
        count = ch_count
        versions = ch_versions
}
