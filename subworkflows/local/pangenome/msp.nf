include { MSPMINER_MSPMINER } from "$projectDir/modules/local/mspminer"
include { MSP_SEQUENCES; MSP_ABUNDANCE } from "$projectDir/modules/local/metagear/utils/post_mspminer"


workflow MSP {

    take:
        gene_catalog
        gene_abundance_count
        gene_abundance_rpkm

    main:

        MSPMINER_MSPMINER ( gene_abundance_count )

        ch_gene_catalog = gene_catalog.map { [ [id: "pangenome"], it[1] ] }
        ch_mspminer_table = MSPMINER_MSPMINER.out.mspminer_main_table.map { [ [id: "pangenome"], it[1] ] }
        ch_post_mspminer = ch_gene_catalog.join(ch_mspminer_table)

        MSP_SEQUENCES ( ch_post_mspminer )

        ch_gene_rpkm = gene_abundance_rpkm.map { [ [id: "pangenome"], it[1] ] }
        ch_msp_abundance = ch_gene_rpkm.join(ch_mspminer_table)

        MSP_ABUNDANCE ( ch_msp_abundance, "median" )

        ch_versions = MSPMINER_MSPMINER.out.versions
                        .mix(MSP_SEQUENCES.out.versions)
                        .mix(MSP_ABUNDANCE.out.versions)

    emit:
        pangenome_dir = MSP_SEQUENCES.out.pangenome_dir
        pangenome_files = MSP_SEQUENCES.out.pangenome_files
        versions = ch_versions
}
