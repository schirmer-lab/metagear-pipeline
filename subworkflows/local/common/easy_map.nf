include { createExistingFileChannel } from "$projectDir/subworkflows/local/utils/existing_data"

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"
include { ABUNDANCE } from "$projectDir/subworkflows/local/common/abundance"

// Map reads against a catalog the caller supplies. No assembly, no gene calling,
// no clustering: the catalog is an input, not something this workflow produces.
// params.catalog_label names the output directory and the merged tables, and
// defaults to the catalog's file name without its extensions.
def catalogLabel( catalog ) {
    if ( params.catalog_label ) return params.catalog_label
    return file(catalog).name.replaceFirst(/\.(fa|fasta|fna|fa\.gz|fasta\.gz|fna\.gz)$/, '')
}

workflow EASY_MAP_INIT {

    main:
        if ( !params.input )   { exit 1, 'Input samplesheet not specified!' }
        if ( !params.catalog ) { exit 1, 'No catalog to map against: set --catalog <fasta>' }

        INPUT_CHECK ( file(params.input), "reads" )

    emit:
        validated_input = INPUT_CHECK.out.validated_input
        versions = INPUT_CHECK.out.versions
}

workflow EASY_MAP {

    take:
        reads // [ meta, reads ]

    main:
        ch_versions = Channel.empty()

        label = catalogLabel( params.catalog )
        ch_catalog = createExistingFileChannel ( params.catalog, { [ [id: label], it ] } )

        ch_abundance_input = reads.combine( ch_catalog )
                                .map { meta_reads, read_files, meta_catalog, catalog ->
                                    [ meta_reads + [label: meta_catalog.id], read_files, catalog ] }

        ABUNDANCE ( ch_abundance_input, 'contig', file("$projectDir/assets/empty.txt") )
        ch_versions = ch_versions.mix( ABUNDANCE.out.versions )

    emit:
        count = ABUNDANCE.out.count
        covered_bases = ABUNDANCE.out.covered_bases
        rpkm = ABUNDANCE.out.rpkm
        tpm = ABUNDANCE.out.tpm
        alignments = ABUNDANCE.out.alignments
        versions = ch_versions
}
