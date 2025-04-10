/* --- IMPORTS --- */
include { METAPHLAN_MAKEDB } from "$projectDir/modules/local/metaphlan4.1/makedb/main"

/* ---  MAIN WORKFLOW --- */
workflow DATABASES {

    main:
        ch_versions = Channel.empty()

        METAPHLAN_MAKEDB ( )

    emit:
        metaphlan_db = METAPHLAN_MAKEDB.out.db
        versions = ch_versions

}
