//
// Summary for workflows (software versions and multiqc)
//

include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from "$projectDir/subworkflows/nf-core/utils_nfcore_pipeline"
include { softwareVersionsToYAML } from "$projectDir/subworkflows/nf-core/utils_nfcore_pipeline"
include { methodsDescriptionText } from "$projectDir/subworkflows/local/utils_nfcore_metagear_pipeline"

/* --- CONFIG FILES --- */
ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo = params.multiqc_logo ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


/* --- IMPORT LOCAL MODULES/SUBWORKFLOWS --- */
include { MULTIQC } from "$projectDir/modules/nf-core/multiqc/main"

workflow SUMMARY {

    take:
        ch_versions
        ch_summary_data

    main:

        //
        // Collate and save software versions
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  ''  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }

        //
        // MODULE: MultiQC
        //
        summary_params = paramsSummaryMap( workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

        ch_methods_description = Channel.value( methodsDescriptionText(ch_multiqc_custom_methods_description) )

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

        ch_multiqc_files = ch_multiqc_files.mix( ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml') )
        ch_multiqc_files = ch_multiqc_files.mix( ch_methods_description.collectFile( name: 'methods_description_mqc.yaml', sort: true ) )
        ch_multiqc_files = ch_multiqc_files.mix( ch_summary_data.collect() )

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )

    emit:
        multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
        versions = ch_versions   // channel: [ path(versions.yml) ]

}
