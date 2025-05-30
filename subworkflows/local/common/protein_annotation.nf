
include { CDHIT_CDHIT } from "$projectDir/modules/local/cdhit/cdhit"
include { SEQKIT_SPLIT2 } from "$projectDir/modules/nf-core/seqkit/split2"
include { INTERPROSCAN } from "$projectDir/modules/local/interproscan/main"
include { FUNCTIONALGROUP_ANNOTATION } from "$projectDir/modules/local/metagear/utils/functional_group_annotation"

workflow PROTEIN_ANNOTATION_INIT {
    main:
        if (params.catalog) {ch_catalog = file(params.catalog)} else { exit 1, 'Input catalog file [fasta format with DNA sequences] not specified!' }

        ch_catalog = Channel.fromPath("${params.catalog}", checkIfExists: true).first()
            .map { it -> [ [id: "gene_catalog"], it] }

    emit:
        catalog_input = ch_catalog
}


workflow PROTEIN_ANNOTATION {

    take:
        input_gene_catalog // [meta, PATH (DNA sequences of gene catalog)]

    main:

        // generate protein catalog
        CDHIT_CDHIT ( ch_cdhit_input )

        // split protein sequences into 1000 fasta files
        SEQKIT_SPLIT2 ( CDHIT_CDHIT.out.fasta )

        // annotate each splited files using interproscan
        INTERPROSCAN ( SEQKIT_SPLIT2.out.reads, "tsv" )

        // create a new channel to collect all interproscan files
        ch_merged_interproscan = INTERPROSCAN.out.tsv.map(it -> it[1]).collect().map(it -> [ ["id": "merged_interproscan_ann"], it])
        // ch_merged_interproscan.view()
        FUNCTIONALGROUP_ANNOTATION ( ch_merged_interproscan )

        ch_versions = CDHIT_CDHIT.out.versions
                        .mix(SEQKIT_SPLIT2.out.versions)
                        .mix(INTERPROSCAN.out.versions.first())
                        .mix(FUNCTIONALGROUP_ANNOTATION.out.versions)

    emit:
        protein_catalog_fasta = CDHIT_CDHIT.out.fasta
        protein_catalog_fasta_clusters = CDHIT_CDHIT.out.clusters
        hits_channel = INTERPROSCAN.out.tsv
        versions = ch_versions
}
