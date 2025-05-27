include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { COVERM_MAKE } from "$projectDir/modules/local/coverm/make"
include { COVERM_CONTIG_BATCH as COVERM_CONTIG; COVERM_CONTIG_MERGE as MERGE_COUNT; COVERM_CONTIG_MERGE as MERGE_RPKM } from "$projectDir/modules/local/coverm/contig"
include { MSPMINER_MSPMINER } from "$projectDir/modules/local/mspminer"


include {BWA_INDEX} from "$projectDir/modules/nf-core/bwa/index"

workflow GENE_PROFILE_INIT {
    main:
        if ( params.input ) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
        if (params.catalog) {ch_catalog = file(params.catalog)} else { exit 1, 'Input catalog file [fasta format] not specified!' }

        INPUT_CHECK ( ch_input, "reads" )

        ch_catalog = Channel.fromPath("${params.catalog}", checkIfExists: true).first()
            .map { it -> [ [id: "gene_catalog"], it] }

    emit:
        validated_input = INPUT_CHECK.out.validated_input
        catalog_input = ch_catalog
        versions = INPUT_CHECK.out.versions
}


workflow GENE_PROFILE {

    take:
        label
        input_reads // [meta, reads]
        assembly_catalog // [ meta, path ]

    main:
        BWA_INDEX ( assembly_catalog ) // build index

        ch_index = BWA_INDEX.out.index.map {[it[1]]}
        ch_make = input_reads.combine(ch_index)

        COVERM_MAKE ( ch_make )

        COVERM_MAKE.out.alignments
            .map { it -> it[1] }
            .collect()
            .map { it -> tuple([id: label], it)}
            .set { ch_coverm_contig }
        
        
        // generate summary table, e.g. count, rpkm, tpm
        // more information about coverm: https://wwood.github.io/CoverM/coverm-contig.html
        COVERM_CONTIG ( ch_coverm_contig )

        ch_coverm_contig_merge_count = COVERM_CONTIG.out.abundance_count.map { it -> tuple( [id: it[0].id + '_count'], it[1] ) }
        MERGE_COUNT ( ch_coverm_contig_merge_count )

        ch_coverm_contig_merge_rpkm = COVERM_CONTIG.out.abundance_rpkm.map { it -> tuple( [id: it[0].id + '_rpkm'], it[1] ) }
        MERGE_RPKM ( ch_coverm_contig_merge_rpkm)

        // run MSPMINER
        MSPMINER_MSPMINER (MERGE_COUNT.out.abundance_count)

        // summary channel version
        ch_versions = BWA_INDEX.out.versions
                        .mix(COVERM_MAKE.out.versions)
                        .mix(COVERM_CONTIG.out.versions)
                        .mix(MSPMINER_MSPMINER.out.versions)

    emit:
        catalog_index = BWA_INDEX.out.index
        alignments = COVERM_MAKE.out.alignments
        abundance_count = MERGE_COUNT.out.abundance_count
        abundance_rpkm = MERGE_RPKM.out.abundance_count
        // abundance_tpm = COVERM_CONTIG.out.abundance_tpm
        // abundance_trimmed_mean = COVERM_CONTIG.out.abundance_trimmed_mean
        mspminer_output_dir = MSPMINER_MSPMINER.out.mspminer_result
        mspminer_main_table = MSPMINER_MSPMINER.out.mspminer_main_table
        versions = ch_versions
}

