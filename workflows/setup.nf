/* --- IMPORTS --- */

include { DATABASES; DATABASES_INIT } from "$projectDir/subworkflows/local/setup/databases"
include { SUMMARY } from "$projectDir/subworkflows/local/common/summary"

/* --- MAIN WORKFLOW --- */
workflow SETUP {

    main:

        ch_versions = Channel.empty()
        summary_data = Channel.empty()

        init = DATABASES_INIT ( ).out
        DATABASES ( init.kneaddata_databases, init.humann_databases, init.database_destinations )

        ch_versions = DATABASES.out.versions

    emit:
        ch_versions

}
