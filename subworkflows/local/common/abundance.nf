include {BWA_INDEX} from "$projectDir/modules/nf-core/bwa/index"

include { COVERM_MAKE } from "$projectDir/modules/local/coverm/make"
include { COVERM_CONTIG_BATCH; COVERM_CONTIG_MERGE } from "$projectDir/modules/local/coverm/contig"

include { COVERM_CONTIG } from "$projectDir/modules/local/coverm/contig"



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

        def chunkCounter = new java.util.concurrent.atomic.AtomicInteger(0)

        COVERM_MAKE.out.alignments
            .map { it[1] }  // keep only BAM path
            .buffer(size: 50, remainder: true ) // emit lists of up to 50 BAMs
            .map { chunk ->
                def idx = chunkCounter.incrementAndGet()
                def partId = "${label}_${String.format('%03d', idx)}"
                tuple([id: partId], chunk)  // => [ [id: label_partNNN], [bam1, ..., bam50] ]
            }
            .set { ch_coverm_contig }

        // generate summary table, e.g. count, rpkm, tpm
        COVERM_CONTIG ( ch_coverm_contig )

        // helper to prepare abundance channels -> [ [id: label_suffix], [files...] ]
        def prepAbundance = { suffix, ch ->
            ch
                .map { it[1] }
                .collect()
                .map { tuple([id: "${label}_${suffix}"], it) }
        }

        ch_coverm_merge = prepAbundance('count', COVERM_CONTIG.out.abundance_count)
                            .concat( prepAbundance('rpkm', COVERM_CONTIG.out.abundance_rpkm) )
                            .concat( prepAbundance('tpm', COVERM_CONTIG.out.abundance_tpm) )

        COVERM_CONTIG_MERGE ( ch_coverm_merge )

        // split merged abundance into separate channels by suffix (avoid AST/into issues)
        ch_tpm   = COVERM_CONTIG_MERGE.out.abundance_merged.filter { it[0].id.endsWith('_tpm') }
        ch_rpkm  = COVERM_CONTIG_MERGE.out.abundance_merged.filter { it[0].id.endsWith('_rpkm') }
        ch_count = COVERM_CONTIG_MERGE.out.abundance_merged.filter { it[0].id.endsWith('_count') }

        // summary channel versions
        ch_versions = BWA_INDEX.out.versions
                        .mix(COVERM_MAKE.out.versions)
                        .mix(COVERM_CONTIG.out.versions)

    emit:
        catalog_index = BWA_INDEX.out.index
        alignments = COVERM_MAKE.out.alignments
        tpm = ch_tpm
        rpkm = ch_rpkm
        count = ch_count
        versions = ch_versions
}

