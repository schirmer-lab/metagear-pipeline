include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"
include { createExistingFileChannel } from "$projectDir/subworkflows/local/utils/existing_data"

include { MSPMINER_MSPMINER } from "$projectDir/modules/local/mspminer"
include { MSP_SEQUENCES; MSP_ABUNDANCE } from "$projectDir/modules/local/metagear/utils/post_mspminer"

include { GTDBTK_CLASSIFYWF } from "$projectDir/modules/local/gtdbtk/classifywf"
include { MSP_METAPHLAN_ANNOTATION } from "$projectDir/modules/local/metagear/utils/msp_metaphlan_annotation"


workflow MSP_INIT {

    main:
        if ( !params.input )                      { exit 1, 'Input samplesheet not specified (--input)' }
        if ( !params.representative_genes )       { exit 1, 'Representative-genes catalog not specified (--representative_genes)' }
        if ( !params.representative_genes_count ) { exit 1, 'Representative-genes counts not specified (--representative_genes_count)' }
        if ( !params.representative_genes_rpkm )  { exit 1, 'Representative-genes RPKM not specified (--representative_genes_rpkm)' }

        ch_input = file(params.input)

        gtdb_tk_db = Channel.fromPath("${params.gtdb_tk_db}", checkIfExists: true)

        representative_genes       = createExistingFileChannel ( params.representative_genes,       { [ [id: "all.genes"],       it ] } )
        representative_genes_count = createExistingFileChannel ( params.representative_genes_count, { [ [id: "all.genes_count"], it ] } )
        representative_genes_rpkm  = createExistingFileChannel ( params.representative_genes_rpkm,  { [ [id: "all.genes_rpkm"],  it ] } )

        metaphlan_db = Channel.empty()
        if ( params.metaphlan_profiles ) {
            metaphlan_profiles = Channel.fromPath("${params.metaphlan_profiles}", checkIfExists: true)
        } else {
            metaphlan_db = Channel.fromPath("${params.metaphlan_db}", checkIfExists: true).first()
            metaphlan_profiles = false
        }

        INPUT_CHECK ( ch_input, "reads" )

    emit:
        validated_input            = INPUT_CHECK.out.validated_input
        representative_genes
        representative_genes_count
        representative_genes_rpkm
        metaphlan_profiles
        metaphlan_db
        gtdb_tk_db
        versions                   = INPUT_CHECK.out.versions
}


workflow MSP {

    take:
        gene_catalog
        gene_abundance_count
        gene_abundance_rpkm
        gtdb_tk_db
        metaphlan_profiles

    main:

        MSPMINER_MSPMINER ( gene_abundance_count )

        ch_gene_catalog = gene_catalog.map { [ [id: "pangenome"], it[1] ] }
        ch_mspminer_table = MSPMINER_MSPMINER.out.mspminer_main_table.map { [ [id: "pangenome"], it[1] ] }
        ch_post_mspminer = ch_gene_catalog.join(ch_mspminer_table)

        MSP_SEQUENCES ( ch_post_mspminer )

        ch_gene_rpkm = gene_abundance_rpkm.map { [ [id: "pangenome"], it[1] ] }
        ch_msp_abundance = ch_gene_rpkm.join(ch_mspminer_table)

        MSP_ABUNDANCE ( ch_msp_abundance, "median" )

        GTDBTK_CLASSIFYWF ( MSP_SEQUENCES.out.pangenome_dir.combine( gtdb_tk_db ), false )

        MSP_METAPHLAN_ANNOTATION ( MSP_ABUNDANCE.out.msp_abundance.combine( metaphlan_profiles ), "v4" )

        ch_versions = MSPMINER_MSPMINER.out.versions
                        .mix(MSP_SEQUENCES.out.versions)
                        .mix(MSP_ABUNDANCE.out.versions)

    emit:
        pangenome_dir = MSP_SEQUENCES.out.pangenome_dir
        pangenome_files = MSP_SEQUENCES.out.pangenome_files
        versions = ch_versions
}
