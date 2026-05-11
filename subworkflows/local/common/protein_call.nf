/* --- Assembly and Gene Calling --- */

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

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
        ch_versions = Channel.empty()

        // translate DNA to protein
        TRANSLATE_DNA2PROT ( gene_catalog )
        ch_versions =  ch_versions.mix(TRANSLATE_DNA2PROT.out.versions)

        ch_proteins = TRANSLATE_DNA2PROT.out.prot_fasta_output

    emit:
        proteins = ch_proteins
        versions = ch_versions
}
