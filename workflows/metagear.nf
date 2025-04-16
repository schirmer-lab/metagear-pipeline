/* --- IMPORT LOCAL SUBWORKFLOWS --- */

include { SETUP } from "$projectDir/workflows/setup"

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { QUALITY_CONTROL_INIT; QUALITY_CONTROL } from "$projectDir/subworkflows/local/common/quality_control"
include { MICROBIAL_PROFILES_INIT; MICROBIAL_PROFILES  } from "$projectDir/subworkflows/local/microbiome/microbial_profiles"


/* --- RUN MAIN WORKFLOW --- */
workflow METAGEAR {

    main:

        ch_versions = Channel.empty()
        ch_summary_data = Channel.empty()

        if (params.workflow == null || params.workflow.trim().isEmpty()) {
            INPUT_CHECK ( file(params.input), "reads" )
        }

        // Setup handler
        if ( params.workflow == "download_databases" ) {
            SETUP ( )
            ch_versions = SETUP.out.versions
        }

        // Quality Control handler
        if ( params.workflow.startsWith("qc_") ) {
            init = QUALITY_CONTROL_INIT ( )
            QUALITY_CONTROL ( init.validated_input, init.kneaddata_refdb )
            ch_versions = QUALITY_CONTROL.out.versions

            ch_summary_data = QUALITY_CONTROL.out.fastqc_zip_pre.collect{it[1]}.ifEmpty([])
                    .mix(QUALITY_CONTROL.out.fastqc_zip_post.collect{it[1]}.ifEmpty([]))
                    .mix(QUALITY_CONTROL.out.summary_plot.collect{it}.ifEmpty([]))
        }

        // Microbial profiles
        if ( params.workflow == "microbial_profiles" ) {
            init = MICROBIAL_PROFILES_INIT ( )
            MICROBIAL_PROFILES ( init.validated_input, init.metaphlan_db, init.uniref90_db, init.chocoplhan_db )
            ch_versions = MICROBIAL_PROFILES.out.versions
        }



    emit:
        versions = ch_versions
        summary_data = ch_summary_data

}
