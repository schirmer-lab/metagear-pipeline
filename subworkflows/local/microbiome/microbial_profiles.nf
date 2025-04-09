/* --- IMPORT LOCAL SUBWORKFLOWS --- */

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"
include { METAPHLAN_METAPHLAN } from "$projectDir/modules/local/metaphlan4.1/metaphlan/main"
include { METAPHLAN_MERGE_PROFILES } from "$projectDir/modules/local/metaphlan4.1/metaphlan/main"
include { SUMMARY } from "$projectDir/subworkflows/local/common/summary"

/* --- RUN MAIN WORKFLOW --- */
workflow MICROBIAL_PROFILES {

    take:
        ch_input // channel: samplesheet read in from --input

    main:

        ch_versions = Channel.empty()
        summary_data = Channel.empty()

        METAPHLAN_METAPHLAN ( ch_input, file(params.metaphlan_db) )

        ch_all_microbial_profiles = METAPHLAN_METAPHLAN.out.microbial_profile
                            .map { [ [id: 'microbial'], it[1] ] }
                            .groupTuple(by: 0)

        METAPHLAN_MERGE_PROFILES( ch_all_microbial_profiles )

        ch_versions = METAPHLAN_METAPHLAN.out.versions.first()
                        .mix( METAPHLAN_MERGE_PROFILES.out.versions.first() )

        SUMMARY ( ch_versions, summary_data )

    emit:
        multiqc_report = SUMMARY.out.multiqc_report

}
