include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { GENE_CALL } from "$projectDir/subworkflows/local/common/gene_call"
include { ABUNDANCE as GENE_ABUNDANCE } from "$projectDir/subworkflows/local/common/abundance"

include { MSPMINER_MSPMINER } from "$projectDir/modules/local/mspminer"
include { TRANSLATE_DNA2PROT } from "$projectDir/modules/local/metagear/utils/translate_dna2prot"

include { PROTEIN_ANNOTATION } from "$projectDir/subworkflows/local/common/protein_annotation"

workflow GENE_ANALYSIS_INIT {

    main:
        if ( params.input ) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

        INPUT_CHECK ( ch_input, "reads" )

    emit:
        validated_input = INPUT_CHECK.out.validated_input
        versions = INPUT_CHECK.out.versions
}


workflow GENE_ANALYSIS {

    take:
        clean_reads // [meta, reads]

    main:

        GENE_CALL ( clean_reads )

        GENE_ABUNDANCE ("label", clean_reads, GENE_CALL.out.genes )

        MSPMINER_MSPMINER ( GENE_ABUNDANCE.out.abundance_count )

        // translate DNA to protein
        TRANSLATE_DNA2PROT ( GENE_CALL.out.genes )
        ch_cdhit_input = TRANSLATE_DNA2PROT.out.prot_fasta_output.map(it -> [[id: "cohort"],it[1]])

        PROTEIN_ANNOTATION ( ch_cdhit_input )

        // summary channel version
        ch_versions = GENE_CALL.out.versions
                        .mix(GENE_ABUNDANCE.out.versions)
                        .mix(MSPMINER_MSPMINER.out.versions)
                        .mix(TRANSLATE_DNA2PROT.out.versions)
                        .mix(PROTEIN_ANNOTATION.out.versions)

    emit:
        // catalog_index = BWA_INDEX.out.index
        // alignments = COVERM_MAKE.out.alignments
        // abundance_count = MERGE_COUNT.out.abundance_count
        // abundance_rpkm = MERGE_RPKM.out.abundance_count
        // abundance_tpm = COVERM_CONTIG.out.abundance_tpm
        // abundance_trimmed_mean = COVERM_CONTIG.out.abundance_trimmed_mean
        // mspminer_output_dir = MSPMINER_MSPMINER.out.mspminer_result
        // mspminer_main_table = MSPMINER_MSPMINER.out.mspminer_main_table
        versions = ch_versions
}
