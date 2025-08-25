/* --- IMPORTS --- */
include { KNEADDATA_DATABASE } from "$projectDir/modules/local/kneaddata/main"
include { METAPHLAN_MAKEDB } from "$projectDir/modules/local/metaphlan4.1/makedb/main"
include { HUMANN_DATABASES } from "$projectDir/modules/local/humann3/main"

include { GENOMAD_DOWNLOAD } from "$projectDir/modules/nf-core/genomad/download/main"
include { CHECKV_DOWNLOADDATABASE } from "$projectDir/modules/nf-core/checkv/downloaddatabase/main"

include { EXPORT_DATABASES } from "$projectDir/modules/local/metagear/export_databases"

/* ---  INITIALIZATION WORKFLOW --- */
workflow DATABASES_INIT {
    main:

        ch_kneaddata_databases = Channel.from( [ ['human_genome', 'bowtie2'] ] )

        ch_humann_databases = Channel.from( ['chocophlan', 'full'], ['uniref', 'uniref90_diamond'] )

        ch_database_destinations = Channel.from( ['metaphlan', file( params.metaphlan_db ) ],
                                            ['chocophlan', file( params.humann3_nucleo ) ],
                                            ['uniref', file( params.humann3_uniref90 ) ],
                                            ['human_genome', file( params.kneaddata_refdb[0] ) ],
                                            ['genomad', file( params.genomad_db ) ],
                                            ['checkv', file( params.checkv_db ) ] )

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

        kneaddata = KNEADDATA_DATABASE( ch_kneaddata_databases )
        ch_versions = ch_versions.mix( kneaddata.versions.first() )

        metaphlan = METAPHLAN_MAKEDB ( )
        ch_versions = ch_versions.mix( metaphlan.versions.first() )

        humann = HUMANN_DATABASES ( ch_humann_databases )
        ch_versions = ch_versions.mix( humann.versions.first() )

        genomad = GENOMAD_DOWNLOAD()
        ch_genomad_database = genomad.genomad_db.map { [ "genomad", it ] }
        ch_versions = ch_versions.mix( genomad.versions )

        checkv = CHECKV_DOWNLOADDATABASE()
        ch_checkv_database = checkv.checkv_db.map { [ "checkv", it ] }
        ch_versions = ch_versions.mix( checkv.versions )

        ch_databases_data = metaphlan.database.concat( humann.database )
                                .concat( kneaddata.database )
                                .concat( ch_genomad_database )
                                .concat( ch_checkv_database )

        ch_databases_data_and_destination = ch_databases_data.join( ch_database_destinations, by: 0 )
                                                .map { [ [id: it[0]], it[1], it[2] ] }


        EXPORT_DATABASES ( ch_databases_data_and_destination )
        ch_versions = ch_versions.mix( EXPORT_DATABASES.out.versions )

    emit:
        versions = ch_versions

}
