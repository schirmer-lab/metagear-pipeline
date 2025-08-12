include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { GENE_CALL } from "$projectDir/subworkflows/local/common/gene_call"
include { PROTEIN_CALL } from "$projectDir/subworkflows/local/common/protein_call"

include { ABUNDANCE as GENE_ABUNDANCE } from "$projectDir/subworkflows/local/common/abundance"

include { PROTEIN_ANNOTATION } from "$projectDir/subworkflows/local/common/protein_annotation"

include { MSP } from "$projectDir/subworkflows/local/pangenome/msp"

include { GTDBTK_CLASSIFYWF } from "$projectDir/modules/local/gtdbtk/classifywf"

workflow GENE_ANALYSIS_INIT {

    main:
        if ( params.input ) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

        gtdb_tk_db = Channel.fromPath("${params.gtdb_tk_db}", checkIfExists: true)

        INPUT_CHECK ( ch_input, "reads" )

    emit:
        validated_input = INPUT_CHECK.out.validated_input
        gtdb_tk_db
        versions = INPUT_CHECK.out.versions
}


workflow GENE_ANALYSIS {

    take:
        clean_reads // [meta, reads]
        gtdb_tk_db

    main:

        GENE_CALL ( clean_reads )

        PROTEIN_CALL ( GENE_CALL.out.gene_catalog )

        GENE_ABUNDANCE ("gene_abundance", clean_reads, GENE_CALL.out.gene_catalog )

        MSP ( GENE_CALL.out.gene_catalog, GENE_ABUNDANCE.out.count, GENE_ABUNDANCE.out.rpkm )

        GTDBTK_CLASSIFYWF ( MSP.out.pangenome_dir.combine( gtdb_tk_db ), false )

        PROTEIN_ANNOTATION ( PROTEIN_CALL.out.protein_catalog )

        // summary channel version
        ch_versions = GENE_CALL.out.versions
                        .mix(PROTEIN_CALL.out.versions)
                        .mix(GENE_ABUNDANCE.out.versions)
                        .mix(MSP.out.versions)
                        .mix(PROTEIN_ANNOTATION.out.versions)


    emit:
        // TODO: implement emission of all relevant channels
        versions = ch_versions
}
