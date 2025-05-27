include { TRANSLATE_DNA2PROT } from "$projectDir/modules/local/metagear/utils/translate_dna2prot"
include { CDHIT_CDHIT } from "$projectDir/modules/local/cdhit/cdhit"
include { SPLIT_FASTA } from "$projectDir/modules/local/metagear/utils/split_fasta"
include { INTERPROSCAN } from "$projectDir/modules/local/interproscan/main"
include { FUNCTIONALGROUP_ANNOTATION } from "$projectDir/modules/local/metagear/utils/functional_group_annotation"

workflow PROTEIN_PROFILE_INIT {
    main:
        if (params.catalog) {ch_catalog = file(params.catalog)} else { exit 1, 'Input catalog file [fasta format with DNA sequences] not specified!' }

        ch_catalog = Channel.fromPath("${params.catalog}", checkIfExists: true).first()
            .map { it -> [ [id: "gene_catalog"], it] }
        // ch_catalog.view()

    emit:
        catalog_input = ch_catalog
}



workflow PROTEIN_PROFILE {

    take:
        input_gene_catalog // [meta, PATH (DNA sequences of gene catalog)]

    main:

        // translate DNA to protein
        TRANSLATE_DNA2PROT ( input_gene_catalog )
        ch_cdhit_input = TRANSLATE_DNA2PROT.out.prot_fasta_output.map(it -> [[id: "cohort"],it[1]])
        // generate protein catalog
        CDHIT_CDHIT ( ch_cdhit_input )

        // split protein sequences into 1000 fasta files
        SPLIT_FASTA ( CDHIT_CDHIT.out.fasta.map(it -> [it[0],it[1],1000]) )
        // 
        SPLIT_FASTA.out.prot_fasta_split
            .flatMap { dir ->
                // Expand all .faa files in that directory
                dir.listFiles().findAll { it.name.endsWith('.faa') }.collect { file ->
                def file_name = file.name.replaceFirst(/\.faa$/, '')
                tuple([id: file_name], file)
                }
            }
            .set { ch_interproscan_input }

        // annotate each splited files using interproscan
        INTERPROSCAN (ch_interproscan_input, "tsv")

        // create a new channel to collect all interproscan files 
        ch_merged_interproscan = INTERPROSCAN.out.tsv.map(it -> it[1]).collect().map(it -> [ ["id": "merged_interproscan_ann"], it])
        ch_merged_interproscan.view()
        FUNCTIONALGROUP_ANNOTATION ( ch_merged_interproscan )

        ch_versions = TRANSLATE_DNA2PROT.out.versions
                        .mix(CDHIT_CDHIT.out.versions)
                        .mix(SPLIT_FASTA.out.versions)
                        .mix(INTERPROSCAN.out.versions.first())
                        .mix(FUNCTIONALGROUP_ANNOTATION.out.versions)

    emit:
        translated_gene_catalog = TRANSLATE_DNA2PROT.out.prot_fasta_output
        protein_catalog_fasta = CDHIT_CDHIT.out.fasta
        protein_catalog_fasta_clusters = CDHIT_CDHIT.out.clusters
        hits_channel = INTERPROSCAN.out.tsv
        versions = ch_versions
}



