/* --- IMPORT LOCAL SUBWORKFLOWS --- */

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"
include { METAPHLAN_METAPHLAN } from "$projectDir/modules/local/metaphlan4.1/metaphlan/main"
include { SUMMARY } from "$projectDir/subworkflows/local/common/summary"

/* --- RUN MAIN WORKFLOW --- */
workflow MICROBIAL_PROFILES {

    take:
        ch_input // channel: samplesheet read in from --input

    main:

        ch_versions = Channel.empty()
        summary_data = Channel.empty()

        // INPUT_CHECK ( ch_input, "reads" )

        METAPHLAN_METAPHLAN ( ch_input, file(params.metaphlan_db) )

        ch_versions = METAPHLAN_METAPHLAN.out.versions.first()

        SUMMARY ( ch_versions, summary_data )

    emit:
        multiqc_report = SUMMARY.out.multiqc_report

}

/* --- END --- */