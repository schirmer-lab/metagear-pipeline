/* --- IMPORTS --- */
include { KNEADDATA_DATABASE } from "$projectDir/modules/local/kneaddata/main"
include { METAPHLAN_MAKEDB } from "$projectDir/modules/local/metaphlan4.1/makedb/main"
include { HUMANN_DATABASES } from "$projectDir/modules/local/humann3/main"

include { EXPORT_DATABASES } from "$projectDir/modules/local/metagear/export_databases"

/* ---  INITIALIZATION WORKFLOW --- */
workflow DATABASES_INIT {
    main:

        ch_kneaddata_databases = Channel.from( [ ['human_genome', 'bowtie2'] ] )

        ch_humann_databases = Channel.from( ['chocophlan', 'full'], ['uniref', 'uniref90_diamond'] )

        ch_database_destinations = Channel.from( ['metaphlan', file( params.metaphlan_db ) ],
                                            ['chocophlan', file( params.humann3_nucleo ) ],
                                            ['uniref', file( params.humann3_uniref90 ) ],
                                            ['human_genome', file( params.kneaddata_refdb[0] ) ] )

        //TODO: Currently only 1 kneaddata database is supported. Ensure ch_kneaddata_databases keep consistent with ch_database_destinations.

    emit:
        kneaddata_databases = ch_kneaddata_databases
        humann_databases = ch_humann_databases
        database_destinations = ch_database_destinations

}


/* ---  MAIN WORKFLOW --- */
workflow DATABASES {
    take:
        ch_kneaddata_databases
        ch_humann_databases
        ch_database_destinations

    main:
        ch_versions = Channel.empty()

        KNEADDATA_DATABASE( ch_kneaddata_databases )

        METAPHLAN_MAKEDB ( )

        HUMANN_DATABASES ( ch_humann_databases )

        ch_databases_data = METAPHLAN_MAKEDB.out.database.concat( HUMANN_DATABASES.out.database )
                                            .concat( KNEADDATA_DATABASE.out.database )

        ch_databases_data_and_destination = ch_databases_data.join( ch_database_destinations, by: 0 )
                                                .map { [ [id: it[0]], it[1], it[2] ] }

        EXPORT_DATABASES ( ch_databases_data_and_destination )

        ch_versions = KNEADDATA_DATABASE.out.versions.first()
                        .mix( METAPHLAN_MAKEDB.out.versions.first() )
                        .mix( HUMANN_DATABASES.out.versions.first() )
                        .mix( EXPORT_DATABASES.out.versions.first() )

    emit:
        versions = ch_versions

}
