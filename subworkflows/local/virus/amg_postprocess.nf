include { EXTRACT_AMGS } from "$projectDir/modules/local/metagear/mge/amg_functions"
include { MMSEQS_EASYSEARCH } from "$projectDir/modules/local/mmseqs/easysearch"
include { MAP_AMG_CATALOG } from "$projectDir/modules/local/metagear/mge/amg_functions"

workflow AMG_POSTPROCESS {

    take:
        amg_info // [meta, amg_summary, amg_fna, amg_faa, tpm, rpkm, count]
        virus_proteins

    main:
        ch_versions = Channel.empty()

        ch_extract_amg = amg_info.map { [ it[0], it[1], it[2], it[3] ] }
        EXTRACT_AMGS ( ch_extract_amg )

        search_ch = EXTRACT_AMGS.out.amgs_faa.combine( virus_proteins )

        MMSEQS_EASYSEARCH ( search_ch )

        amg_catalog_input = amg_info
                            .map { [ it[0], it[1], it[4], it[5], it[6] ] }
                            .join( MMSEQS_EASYSEARCH.out.tsv )


        MAP_AMG_CATALOG ( amg_catalog_input )

    emit:
        amgs_faa = EXTRACT_AMGS.out.amgs_faa
        amgs_fna = EXTRACT_AMGS.out.amgs_fna
        versions = ch_versions
}
