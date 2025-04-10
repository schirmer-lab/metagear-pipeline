/* --- IMPORTS --- */

include { DATABASES } from "$projectDir/subworkflows/local/setup/databases"
include { SUMMARY } from "$projectDir/subworkflows/local/common/summary"

/* --- MAIN WORKFLOW --- */
workflow SETUP {

    main:

        ch_versions = Channel.empty()
        summary_data = Channel.empty()

        DATABASES ( )

        SUMMARY ( ch_versions, summary_data )

    emit:
        multiqc_report = SUMMARY.out.multiqc_report

}
