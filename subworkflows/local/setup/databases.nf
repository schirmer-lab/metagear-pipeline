/* --- IMPORTS --- */
include { METAPHLAN_MAKEDB } from "$projectDir/modules/local/metaphlan4.1/makedb/main"
include { HUMANN_DATABASES } from "$projectDir/modules/local/humann3/main"
include { EXPORT_DATABASES } from "$projectDir/modules/local/metagear/export_databases"

/* ---  MAIN WORKFLOW --- */
workflow DATABASES {

    main:
        ch_versions = Channel.empty()

        METAPHLAN_MAKEDB ( )

        ch_humann_databases = Channel.from( ['chocophlan', 'full'], ['uniref', 'uniref90_diamond'] )

        HUMANN_DATABASES ( ch_humann_databases )

        // EXPORT_DATABASES()
        ch_databases_locations = Channel.from( ['metaphlan', file( params.metaphlan_db ) ],
                                            ['chocophlan', file( params.humann3_nucleo ) ],
                                            ['uniref', file( params.humann3_uniref90 ) ] )

        ch_databases_data = METAPHLAN_MAKEDB.out.database.concat( HUMANN_DATABASES.out.database )

        ch_databases_export = ch_databases_data.join( ch_databases_locations, by: 0 )
                                                .map { [ [id: it[0]], it[1], it[2] ] }

        EXPORT_DATABASES ( ch_databases_export )

    emit:
        metaphlan_db = METAPHLAN_MAKEDB.out.database
        versions = ch_versions

}
