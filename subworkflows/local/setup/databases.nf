/* --- IMPORTS --- */
include { KNEADDATA_DATABASE } from "$projectDir/modules/local/kneaddata/main"
include { METAPHLAN_MAKEDB } from "$projectDir/modules/local/metaphlan4.1/makedb/main"
include { HUMANN_DATABASES } from "$projectDir/modules/local/humann3/main"
include { GTDBTK_DOWNLOAD_DB } from "$projectDir/modules/local/gtdbtk/download/main"

include { GENOMAD_DOWNLOAD } from "$projectDir/modules/nf-core/genomad/download/main"
include { CHECKV_DOWNLOADDATABASE } from "$projectDir/modules/nf-core/checkv/downloaddatabase/main"
include { CHECKM2_DATABASEDOWNLOAD } from "$projectDir/modules/nf-core/checkm2/databasedownload/main"
include { MMSEQS_DATABASES } from "$projectDir/modules/local/mmseqs/databases/main"
// include { VIRSORTER2_SETUP } from "$projectDir/modules/local/virsorter2/setup"
include { DRAM_SETUP } from "$projectDir/modules/local/dram/setup"
include { IPHOP_DOWNLOAD } from "$projectDir/modules/local/iphop/download/main"
include { PHAROKKA_INSTALLDATABASES } from "$projectDir/modules/nf-core/pharokka/installdatabases/main"

include { AMRFINDERPLUS_UPDATE } from "$projectDir/modules/nf-core/amrfinderplus/update/main"

include { PHOLD_INSTALL } from "$projectDir/modules/local/phold/install/main"

include { EXPORT_DATABASES } from "$projectDir/modules/local/metagear/export_databases"

/* ---  INITIALIZATION WORKFLOW --- */
workflow DATABASES_INIT {
    main:

        ch_kneaddata_databases = Channel.from( [ ['human_genome', 'bowtie2'] ] )

        ch_humann_databases = Channel.from( ['chocophlan', 'full'], ['uniref', 'uniref90_diamond'] )

        ch_database_destinations = Channel.from( [ 'metaphlan', file( params.metaphlan_db ) ],
                                            [ 'chocophlan', file( params.humann3_nucleo ) ],
                                            [ 'uniref', file( params.humann3_uniref90 ) ],
                                            [ 'human_genome', file( params.kneaddata_refdb[0] ) ],
                                            [ 'gtdb_tk', file( params.gtdb_tk_db ) ],
                                            [ 'genomad', file( params.genomad_db ) ],
                                            [ 'checkv', file( params.checkv_db ) ],
                                            // [ 'virsorter2', file( params.virsorter2_db ) ],
                                            [ 'dram', file( params.dram_db ) ],
                                            [ 'iphop', file( params.iphop_db ) ],
                                            [ 'amrfinder', file( params.amrfinder_db ) ],
                                            [ 'pharokka', file( params.pharokka_db ) ],
                                            [ 'checkm2', file( params.checkm2_db ) ],
                                            [ 'mmseqs_taxonomy', file( params.mmseqs_taxonomy_db ) ],
                                            [ 'phold', file( params.phold_db ?: '/dev/null' ) ] )

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
        ch_databases_data = Channel.empty()

        if ( params.databases == "all" || params.databases.contains("genes") ) {

            kneaddata = KNEADDATA_DATABASE( ch_kneaddata_databases )
            ch_versions = ch_versions.mix( kneaddata.versions )

            metaphlan = METAPHLAN_MAKEDB ( )
            ch_versions = ch_versions.mix( metaphlan.versions )

            humann = HUMANN_DATABASES ( ch_humann_databases )
            ch_versions = ch_versions.mix( humann.versions )

            gtdbtk = GTDBTK_DOWNLOAD_DB ( )

            gtdbtk_database = gtdbtk.database.map { [ "gtdb_tk", it ] }
            ch_versions = ch_versions.mix( gtdbtk.versions )

            ch_databases_data = ch_databases_data
                                .concat( kneaddata.database )
                                .concat( humann.database )
                                .concat( metaphlan.database )
                                .concat( gtdbtk.database )
        }

        if ( params.databases == "all" || params.databases.contains("virus") ) {

            genomad = GENOMAD_DOWNLOAD ( )
            ch_genomad_database = genomad.genomad_db.map { [ "genomad", it ] }
            ch_versions = ch_versions.mix( genomad.versions )
            ch_databases_data = ch_databases_data.concat( ch_genomad_database )

            checkv = CHECKV_DOWNLOADDATABASE ( )
            ch_checkv_database = checkv.checkv_db.map { [ "checkv", it ] }
            ch_versions = ch_versions.mix( checkv.versions )
            ch_databases_data = ch_databases_data.concat( ch_checkv_database )

            // virsorter2 = VIRSORTER2_SETUP ( )
            // virsorter2_database = virsorter2.virsorter2_db.map { [ "virsorter2", it ] }
            // ch_versions = ch_versions.mix( virsorter2.versions )

            dram = DRAM_SETUP ( )
            dram_database = dram.dram_db.map { [ "dram", it ] }
            // ch_versions = ch_versions.mix( dram.versions ) //TODO: Fix version.yml from DRAM_SETUP
            ch_databases_data = ch_databases_data.concat( dram_database )

            iphop = IPHOP_DOWNLOAD ( )
            iphop_database = iphop.iphop_db.map { [ "iphop", it ] }
            ch_versions = ch_versions.mix( iphop.versions )
            ch_databases_data = ch_databases_data.concat( iphop_database )

            amrfinderplus = AMRFINDERPLUS_UPDATE ( )
            amrfinderplus_database = amrfinderplus.db.map { [ "amrfinder", it ] }
            ch_versions = ch_versions.mix( amrfinderplus.versions )
            ch_databases_data = ch_databases_data.concat( amrfinderplus_database )

            pharokka = PHAROKKA_INSTALLDATABASES ( )
            pharokka_database = pharokka.pharokka_db.map { [ "pharokka", it ] }
            ch_versions = ch_versions.mix( pharokka.versions )
            ch_databases_data = ch_databases_data.concat( pharokka_database )

        }

        if ( params.databases == "all" || params.databases.contains("contig_classification") ) {

            checkm2 = CHECKM2_DATABASEDOWNLOAD ( Channel.value([]) )
            ch_checkm2_database = checkm2.database.map { meta, file -> [ "checkm2", file ] }
            // TODO: CHECKM2_DATABASEDOWNLOAD reports versions via the new topic channel
            //       pattern (not path 'versions.yml'); wire into ch_versions once the
            //       pipeline standardises on topic channels.
            ch_databases_data = ch_databases_data.concat( ch_checkm2_database )

        }

        // PHOLD structural DB (consumed by the `structures` workflow).
        // Standalone group so users can install it without re-checking the
        // viral DB stack; included in "all" for greenfield installs. ~7.7 GB.
        if ( params.databases == "all" || params.databases.contains("structures") ) {

            phold = PHOLD_INSTALL ( )
            ch_phold_database = phold.phold_db.map { [ "phold", it ] }
            ch_versions = ch_versions.mix( phold.versions )
            ch_databases_data = ch_databases_data.concat( ch_phold_database )

        }

        // MMseqs2 taxonomy DB (used by CLASSIFICATION's contig fallback).
        // Triggered both by the umbrella "contig_classification" group AND by the
        // dedicated "mmseqs_taxonomy" group, so `--databases mmseqs_taxonomy` pulls
        // only this one DB.
        if ( params.databases == "all"
                || params.databases.contains("contig_classification")
                || params.databases.contains("mmseqs_taxonomy") ) {

            // Default to GTDB; override with `process.withName:'.*MMSEQS_DATABASES'.ext.db_name`
            // or by setting params.mmseqs_taxonomy_source.
            def mmseqs_source = params.mmseqs_taxonomy_source ?: 'GTDB'
            mmseqs_tax = MMSEQS_DATABASES ( Channel.value(mmseqs_source) )
            ch_mmseqs_taxonomy_database = mmseqs_tax.database.map { dir -> [ "mmseqs_taxonomy", dir ] }
            ch_versions = ch_versions.mix( mmseqs_tax.versions )
            ch_databases_data = ch_databases_data.concat( ch_mmseqs_taxonomy_database )

        }

        ch_databases_data_and_destination = ch_databases_data.join( ch_database_destinations, by: 0 )
                                                .map { [ [id: it[0]], it[1], it[2] ] }

        // ch_databases_data_and_destination.view()

        EXPORT_DATABASES ( ch_databases_data_and_destination )
        // ch_versions = ch_versions.mix( EXPORT_DATABASES.out.versions )

    emit:
        versions = ch_versions

}
