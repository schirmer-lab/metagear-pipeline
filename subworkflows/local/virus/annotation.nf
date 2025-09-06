include { VIRSORTER2 as VIRSORTER2_4DRAMV } from "$projectDir/modules/local/virsorter2/main"
include { DRAMV } from "$projectDir/modules/local/dram/main"
include { IPHOP_PREDICT } from "$projectDir/modules/local/iphop/predict/main"
// include { AMRFINDERPLUS_RUN } from "$projectDir/modules/nf-core/amrfinderplus/run/main"
include { SEQKIT_SPLIT2 } from "$projectDir/modules/nf-core/seqkit/split2"

workflow VIRAL_ANNOTATION_INIT {

    main:
        if (params.input) {ch_catalog = file(params.input)} else { exit 1, 'Input catalog file [fasta format with DNA sequences] not specified!' }

        ch_catalog = Channel.fromPath("${params.input}", checkIfExists: true)
            .map { it -> [ [id: "catalog"], it] }

        ch_protein_catalog = Channel.fromPath("${params.protein_input}", checkIfExists: true)
            .map { it -> [ [id: "protein_catalog", is_proteins: true], it] }

        virsorter2_db = Channel.fromPath("${params.virsorter2_db}", checkIfExists: true)
        dram_db = Channel.fromPath("${params.dram_db}", checkIfExists: true)
        iphop_db = Channel.fromPath("${params.iphop_db}", checkIfExists: true)

    emit:
        catalog_input = ch_catalog
        protein_catalog = ch_protein_catalog
        virsorter2_db
        dram_db
        iphop_db
}


workflow VIRAL_ANNOTATION {

    take:
        contigs // [meta, fasta]
        virsorter2_db
        dram_db
        iphop_db

    main:
        ch_versions = Channel.empty()

        ch_split = contigs.map { meta, path ->
            def newMeta = meta.clone()
            newMeta.single_end = true
            return tuple(newMeta, path)
        }

        // split protein sequences into 1000 fasta files
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
            .set { ch_virsorter2_chunks }

        VIRSORTER2_4DRAMV( ch_virsorter2_chunks, virsorter2_db )
        ch_versions = ch_versions.mix(VIRSORTER2_4DRAMV.out.versions)


        ch_dramv_input = VIRSORTER2_4DRAMV.out.vs2_4dram_virus
                            .join ( VIRSORTER2_4DRAMV.out.vs2_4dra_affi, by: 0)

        DRAMV ( ch_dramv_input, dram_db )
        ch_versions = ch_versions.mix(DRAMV.out.versions)

        IPHOP_PREDICT ( ch_virsorter2_chunks, iphop_db )

    emit:

        amgs = DRAMV.out.amg_summary
        iphop_genus = IPHOP_PREDICT.out.iphop_genus
        iphop_genomes = IPHOP_PREDICT.out.iphop_genome
        versions = ch_versions
}

