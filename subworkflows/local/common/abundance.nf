include { BWA_INDEX } from "$projectDir/modules/nf-core/bwa/index"
include { BWA_INDEX_CHECK } from "$projectDir/modules/local/bwa/check_index"

include { COVERM_MAKE } from "$projectDir/modules/local/coverm/make"
include { COVERM_CONTIG; COVERM_CONTIG_MERGE } from "$projectDir/modules/local/coverm/contig"
include { COVERM_GENOME } from "$projectDir/modules/local/coverm/genome"


workflow ABUNDANCE {

    take:
        reads_with_sequences // [ meta, reads, catalog ]: meta must contain id and label
        mode                 // 'contig' | 'genome' — selects coverm subcommand
        genome_definition    // path: contig→genome TSV consumed when mode == 'genome'.
                             // In contig mode pass a placeholder (e.g.
                             // file("$projectDir/assets/empty.txt")); it is never read.
        existing_index       // directory holding a bwa index for the catalog, or null to
                             // build one. Built for a single catalog: with several, the
                             // check below fails for the ones it does not match.

    main:

        /* -- Build index for abundance estimation (only once) --- */
        ch_catalogs = reads_with_sequences
                        .map { meta, reads, catalog -> [ [id: meta.label, src: catalog.getName()], catalog ] }
                        .unique()

        if ( existing_index ) {
            BWA_INDEX_CHECK ( ch_catalogs.map { meta, catalog ->
                                  [ meta, file(existing_index, checkIfExists: true), catalog ] } )
            ch_index = BWA_INDEX_CHECK.out.index
            ch_index_versions = Channel.empty()
        } else {
            BWA_INDEX ( ch_catalogs )
            ch_index = BWA_INDEX.out.index
            ch_index_versions = BWA_INDEX.out.versions
        }

        // Combine back the index with reads and catalog for abundance estimation
        ch_reads_with_index = reads_with_sequences
                        .map { meta, reads, catalog -> [ [id: meta.label, src: catalog.getName()], meta, reads ] }
                        .combine( ch_index, by: 0 )
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

        if ( mode == 'genome' ) {
            // Attach the cohort-global contig→genome definition to every batch.
            // `combine` against a value/queue channel of one item broadcasts.
            COVERM_GENOME ( ch_coverm_contig.combine( genome_definition ) )
            ch_abundance = COVERM_GENOME.out.abundance
            ch_abund_versions = COVERM_GENOME.out.versions
        } else {
            COVERM_CONTIG ( ch_coverm_contig )
            ch_abundance = COVERM_CONTIG.out.abundance
            ch_abund_versions = COVERM_CONTIG.out.versions
        }

        // One merge per label and metric, the metric taken from the file name.
        ch_coverm_merge = ch_abundance
                            .transpose()
                            .map { meta, tsv ->
                                def metric = ( tsv.name =~ /\.abundance_(.+)\.tsv$/ )[0][1]
                                [ [id: "${meta.label}_${metric}"], tsv ]
                            }
                            .groupTuple (by: 0)

        COVERM_CONTIG_MERGE ( ch_coverm_merge )

        def merged_metric = { metric ->
            COVERM_CONTIG_MERGE.out.abundance_merged.filter { it[0].id.endsWith("_${metric}") }
        }

        // summary channel versions
        ch_versions = COVERM_MAKE.out.versions
                        .mix(ch_abund_versions)
                        .mix(ch_index_versions)
        // ch_versions = Channel.empty()

    emit:
        index = ch_index
        alignments = COVERM_MAKE.out.alignments
        tpm = merged_metric('tpm')
        rpkm = merged_metric('rpkm')
        count = merged_metric('count')
        covered_bases = merged_metric('covered_bases')
        versions = ch_versions
}
