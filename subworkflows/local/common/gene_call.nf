/* --- Assembly and Gene Calling --- */

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { PRODIGAL } from "$projectDir/modules/nf-core/prodigal"
include { FILTER_PRODIGAL } from "$projectDir/modules/local/metagear/utils/filter_prodigal"

include { VAMB_CONCATENATE_FASTA } from "$projectDir/modules/local/vamb/main"

include { CDHIT_CDHITEST } from "$projectDir/modules/local/cdhit/cdhitest/main"

/* --- Initialization for standalone process --- */
workflow GENE_CALL_INIT {
    main:

        if ( params.input ) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

        INPUT_CHECK ( ch_input, "reads" )

    emit:
        validated_input = INPUT_CHECK.out.validated_input
        versions = INPUT_CHECK.out.versions

}

/* --- Main Workflow --- */
workflow GENE_CALL {

    take:
        sequences // meta, reads

    main:

        PRODIGAL ( sequences, "gff" )

        FILTER_PRODIGAL ( PRODIGAL.out.nucleotide_fasta )

        ch_merged_genes = FILTER_PRODIGAL.out.filtered_fasta.map{ it -> it[1] }
                .collect()
                .map{ it -> [ [id: "merged_genes"], it ] }

        VAMB_CONCATENATE_FASTA ( ch_merged_genes )

        CDHIT_CDHITEST ( VAMB_CONCATENATE_FASTA.out.catalog )

        ch_versions = PRODIGAL.out.versions.first()
                        .mix(FILTER_PRODIGAL.out.versions.first())
                        .mix(VAMB_CONCATENATE_FASTA.out.versions)
                        .mix(CDHIT_CDHITEST.out.versions)


    emit:
        genes = CDHIT_CDHITEST.out.fasta
        versions = ch_versions
}
