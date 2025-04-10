/* --- IMPORT LOCAL SUBWORKFLOWS --- */

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { METAPHLAN_METAPHLAN } from "$projectDir/modules/local/metaphlan4.1/metaphlan/main"
include { METAPHLAN_MERGE_PROFILES } from "$projectDir/modules/local/metaphlan4.1/metaphlan/main"

include { HUMANN_FUNCTION; HUMANN_MERGE_PROFILES } from "$projectDir/modules/local/humann3/main"


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

        ch_reads_profiles = ch_input.join (METAPHLAN_METAPHLAN.out.microbial_profile, by: 0)

        HUMANN_FUNCTION ( ch_reads_profiles, file(params.humann3_uniref90), file(params.humann3_nucleo) )

        ch_all_gene_families = HUMANN_FUNCTION.out.gene_family
                                .map { [ [id: 'gene_families'], it[1] ] }
                                .groupTuple(by: 0)

        ch_all_path_abundances = HUMANN_FUNCTION.out.path_abundance
                                .map { [ [id: 'path_abundances'], it[1] ] }
                                .groupTuple(by: 0)

        HUMANN_MERGE_PROFILES ( ch_all_gene_families.concat( ch_all_path_abundances ) )

        ch_versions = METAPHLAN_METAPHLAN.out.versions.first()
                        .mix( METAPHLAN_MERGE_PROFILES.out.versions.first() )
                        .mix( HUMANN_FUNCTION.out.versions.first() )

        SUMMARY ( ch_versions, summary_data )

    emit:
        multiqc_report = SUMMARY.out.multiqc_report

}
