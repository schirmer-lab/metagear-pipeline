include { SEQKIT_SPLIT2 } from "$projectDir/modules/nf-core/seqkit/split2"
include { INTERPROSCAN } from "$projectDir/modules/local/interproscan/main"
include { FUNCTIONALGROUP_ANNOTATION } from "$projectDir/modules/local/metagear/utils/functional_group_annotation"

workflow PROTEIN_ANNOTATION_INIT {
    main:
        if (params.protein_catalog) {ch_catalog = file(params.protein_catalog)} else { exit 1, 'Input catalog file [fasta format with DNA sequences] not specified!' }

        ch_catalog = Channel.fromPath("${params.protein_catalog}", checkIfExists: true).first()
            .map { it -> [ [id: "gene_catalog"], it] }

    emit:
        catalog_input = ch_catalog
}


workflow PROTEIN_ANNOTATION {

    take:
        protein_catalog // [meta, PATH (DNA sequences of gene catalog)]

    main:

        ch_split = protein_catalog.map { meta, path ->
            def newMeta = meta.clone()
            newMeta.single_end = true
            return tuple(newMeta, path)
        }

        // split protein sequences into chunks with n sequences (see config, default 5K)
        SEQKIT_SPLIT2 ( ch_split )

        SEQKIT_SPLIT2.out.reads
            .flatMap { meta, gz ->
                def files = (gz instanceof java.nio.file.Path) ? [gz] : (gz as List)
                files.collect { f ->
                    def fn = f.getFileName().toString()
                    def chunkId = fn.replaceFirst(/\.faa\.gz$/, '')
                    tuple([ id: "${chunkId}" ], f)   // keep full path as Path
                }
            }
            .set { ch_interproscan_input }

        INTERPROSCAN ( ch_interproscan_input, "tsv" )

        // create a new channel to collect all interproscan files
        ch_merged_interproscan = INTERPROSCAN.out.tsv.map( it -> it[1] ).collect().map(it -> [ ["id": "protein_catalog"], it])

        FUNCTIONALGROUP_ANNOTATION ( ch_merged_interproscan )

        ch_versions = SEQKIT_SPLIT2.out.versions
                        .mix(INTERPROSCAN.out.versions.first())
                        .mix(FUNCTIONALGROUP_ANNOTATION.out.versions)

    emit:
        hits_channel = INTERPROSCAN.out.tsv
        versions = ch_versions
}
