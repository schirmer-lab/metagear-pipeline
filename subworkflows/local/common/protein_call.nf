/* --- Assembly and Gene Calling --- */

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

// include { CDHIT_CDHIT } from "$projectDir/modules/local/cdhit/cdhit"
include { MMSEQS_EASY_CLUSTER } from "$projectDir/modules/local/mmseqs/easy_cluster/main"

include { TRANSLATE_DNA2PROT } from "$projectDir/modules/local/metagear/utils/translate_dna2prot"

/* --- Initialization for standalone process --- */
workflow PROTEIN_CALL_INIT {
    main:
        if (params.gene_catalog) {ch_catalog = file(params.gene_catalog)} else { exit 1, 'Input catalog file [fasta format with DNA sequences] not specified!' }

        ch_catalog = Channel.fromPath("${params.gene_catalog}", checkIfExists: true).first()
            .map { it -> [ [id: "gene_catalog"], it] }

    emit:
        gene_catalog = ch_catalog
}


/* --- Main Workflow --- */
workflow PROTEIN_CALL {

    take:
        gene_catalog // meta, sequences

    main:

        // translate DNA to protein
        TRANSLATE_DNA2PROT ( gene_catalog )

        ch_protein_catalog_input = TRANSLATE_DNA2PROT.out.prot_fasta_output.map(it -> [ [id: it[0].id.replace(".genes", ".proteins") ], it[1]])

        // CDHIT_CDHIT ( ch_protein_catalog_input )
        MMSEQS_EASY_CLUSTER ( ch_protein_catalog_input )

        ch_protein_catalog = MMSEQS_EASY_CLUSTER.out.representatives
        ch_protein_catalog_clusters = MMSEQS_EASY_CLUSTER.out.clusters_tsv

        ch_versions = TRANSLATE_DNA2PROT.out.versions.first()
                        .mix(MMSEQS_EASY_CLUSTER.out.versions.first())


    emit:
        protein_catalog = ch_protein_catalog
        protein_catalog_clusters = ch_protein_catalog_clusters
        versions = ch_versions
}
