/* --- IMPORT LOCAL SUBWORKFLOWS --- */

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { METAPHLAN_METAPHLAN } from "$projectDir/modules/local/metaphlan4.1/metaphlan/main"
include { METAPHLAN_MERGE_PROFILES } from "$projectDir/modules/local/metaphlan4.1/metaphlan/main"

include { HUMANN_FUNCTION; HUMANN_MERGE_PROFILES } from "$projectDir/modules/local/humann3/main"


include { SUMMARY } from "$projectDir/subworkflows/local/common/summary"

/* --- INITIALIZATION FOR STANDALONE RUN --- */
workflow MICROBIAL_PROFILES_INIT {

    main:

        if ( params.input ) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

        INPUT_CHECK ( file(params.input), "reads" )
        ch_versions = INPUT_CHECK.out.versions.first()

        metaphlan_db = Channel.fromPath("${params.metaphlan_db}", checkIfExists: true).first()
        uniref90_db = Channel.fromPath("${params.humann3_uniref90}", checkIfExists: true).first()
        chocoplhan_db = Channel.fromPath("${params.humann3_nucleo}", checkIfExists: true).first()

    emit:
        validated_input = INPUT_CHECK.out.validated_input
        metaphlan_db
        uniref90_db
        chocoplhan_db
        versions = ch_versions
}

/* --- MAIN WORKFLOW --- */
workflow MICROBIAL_PROFILES {

    take:
        validated_input // channel: validated_input from INPUT_CHECK or Upstream workflow
        metaphlan_db
        humann3_uniref90_db
        humann3_chocoplhan_db

    main:

        ch_versions = Channel.empty()

        METAPHLAN_METAPHLAN ( validated_input, metaphlan_db )

        ch_all_microbial_profiles = METAPHLAN_METAPHLAN.out.microbial_profile
                                    .map { [ [id: 'microbial'], it[1] ] }
                                    .groupTuple(by: 0)

        METAPHLAN_MERGE_PROFILES( ch_all_microbial_profiles )

        ch_reads_profiles = validated_input.join (METAPHLAN_METAPHLAN.out.microbial_profile, by: 0)

        HUMANN_FUNCTION ( ch_reads_profiles, humann3_uniref90_db, humann3_chocoplhan_db )

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

    emit:
        versions = ch_versions

}
