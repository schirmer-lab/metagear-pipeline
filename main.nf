#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    schirmer-lab/metagear
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/schirmer-lab/metagear
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { METAGEAR  } from "$projectDir/workflows/metagear"
include { SETUP  } from "$projectDir/workflows/setup"
include { PIPELINE_INITIALISATION } from "$projectDir/subworkflows/local/utils_nfcore_metagear_pipeline"
include { PIPELINE_COMPLETION } from "$projectDir/subworkflows/local/utils_nfcore_metagear_pipeline"

include { MICROBIAL_PROFILES } from "$projectDir/subworkflows/local/microbiome/microbial_profiles"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// WORKFLOW: Run main analysis pipeline depending on type of input
workflow SCHIRMERLAB_METAGEAR {

    take:
        samplesheet // channel: samplesheet read in from --input

    main:

        ch_multiqc_report = Channel.empty()

        if ( params.module != "" ) {

        }
        else if ( params.subworkflow != "" ) {

            if ( params.subworkflow == "microbial_profiles" ) {
                MICROBIAL_PROFILES ( samplesheet )
                ch_multiqc_report = MICROBIAL_PROFILES.out.multiqc_report
            }

        }
        else if( params.workflow != "" ) {
            if ( params.workflow == "setup" ) {
                SETUP ( )
                // ch_multiqc_report = MICROBIAL_PROFILES.out.multiqc_report
            }
        }
        else{
            // WORKFLOW: Run pipeline
            METAGEAR ( samplesheet )
            ch_multiqc_report = METAGEAR.out.multiqc_report
        }



    emit:
        multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
        //
        // SUBWORKFLOW: Run initialisation tasks
        //
        PIPELINE_INITIALISATION (
            params.version,
            params.validate_params,
            params.monochrome_logs,
            args,
            params.outdir,
            params.input
        )

        //
        // WORKFLOW: Run main workflow
        //
        SCHIRMERLAB_METAGEAR (
            PIPELINE_INITIALISATION.out.samplesheet
        )
        //
        // SUBWORKFLOW: Run completion tasks
        //
        PIPELINE_COMPLETION (
            params.email,
            params.email_on_fail,
            params.plaintext_email,
            params.outdir,
            params.monochrome_logs,
            params.hook_url,
            SCHIRMERLAB_METAGEAR.out.multiqc_report
        )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
