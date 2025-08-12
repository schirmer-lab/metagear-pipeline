/* --- Assembly and Gene Calling --- */

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { CDHIT_CDHIT } from "$projectDir/modules/local/cdhit/cdhit"
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

        ch_protein_catalog = TRANSLATE_DNA2PROT.out.prot_fasta_output.map(it -> [[id: "protein_catalog"],it[1]])

        CDHIT_CDHIT ( ch_protein_catalog )

        ch_versions = TRANSLATE_DNA2PROT.out.versions.first()
                        .mix(CDHIT_CDHIT.out.versions.first())


    emit:
        protein_catalog = CDHIT_CDHIT.out.fasta
        protein_catalog_clusters = CDHIT_CDHIT.out.clusters
        versions = ch_versions
}
