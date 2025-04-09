/* --- IMPORT LOCAL SUBWORKFLOWS --- */

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"
include { SUMMARY } from "$projectDir/subworkflows/local/common/summary"

/* --- RUN MAIN WORKFLOW --- */
workflow METAGEAR {

    take:
        ch_input // channel: samplesheet read in from --input

    main:

        ch_versions = Channel.empty()
        summary_data = Channel.empty()

        INPUT_CHECK ( ch_input, "reads" )

        ch_versions = ch_versions.mix( INPUT_CHECK.out.versions.first() )
        SUMMARY ( ch_versions, summary_data )

    emit:
        multiqc_report = SUMMARY.out.multiqc_report

}
