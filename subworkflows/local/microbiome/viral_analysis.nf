include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { createExistingDirChannel; createExistingFileChannel } from "$projectDir/subworkflows/local/utils/existing_data"

include { ASSEMBLY } from "$projectDir/subworkflows/local/common/assembly"

include { VIRAL_DETECTION } from "$projectDir/subworkflows/local/virus/detection"
include { VIRAL_ANNOTATION } from "$projectDir/subworkflows/local/virus/annotation"

include {BWA_INDEX} from "$projectDir/modules/nf-core/bwa/index"

include { GENE_CALL; VIRAL_GENE_CALL } from "$projectDir/subworkflows/local/common/gene_call"
include { PROTEIN_CALL } from "$projectDir/subworkflows/local/common/protein_call"
include { PROTEIN_ANNOTATION } from "$projectDir/subworkflows/local/common/protein_annotation"

include { ABUNDANCE } from "$projectDir/subworkflows/local/common/abundance"

include { CLUSTER_SEQUENCES } from "$projectDir/subworkflows/local/common/clustering"

include { COLLECT_TABLES } from "$projectDir/modules/local/metagear/mge/summarize"

workflow VIRAL_ANALYSIS_INIT {

    main:
        if ( params.input ) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

        genomad_db = Channel.fromPath("${params.genomad_db}", checkIfExists: true).first()
        checkv_db = Channel.fromPath("${params.checkv_db}", checkIfExists: true).first()
        virsorter2_db = Channel.fromPath("${params.virsorter2_db}", checkIfExists: true).first()
        dram_db = Channel.fromPath("${params.dram_db}", checkIfExists: true).first()
        iphop_db = Channel.fromPath("${params.iphop_db}", checkIfExists: true).first()
        amrfinder_db = Channel.fromPath("${params.amrfinder_db}", checkIfExists: true).first()

        INPUT_CHECK ( ch_input, "reads" )

    emit:
        genomad_db
        checkv_db
        virsorter2_db
        dram_db
        iphop_db
        amrfinder_db

        reads = INPUT_CHECK.out.validated_input
        versions = INPUT_CHECK.out.versions
}


workflow VIRAL_ANALYSIS {

    take:
        reads // [meta, reads(fast1, fas2)]
        genomad_db
        checkv_db
        virsorter2_db
        dram_db
        iphop_db
        amrfinder_db

    main:

        ch_versions = Channel.empty()

        ch_contigs = Channel.empty()

        if ( params.contigs_dir ) {
            ch_contigs = createExistingDirChannel ( params.contigs_dir, "*.contigs.fa.gz", ".contigs.fa", false )

        }else {
            // Assemble reads into contigs
            ASSEMBLY ( reads )
            ch_versions = ch_versions.mix( ASSEMBLY.out.versions.first() )

            ch_contigs = ASSEMBLY.out.contigs
        }

        VIRAL_DETECTION ( ch_contigs, genomad_db, checkv_db )

        ch_genes = Channel.empty()
        if ( params.genes_dir ) {
            ch_genes = createExistingDirChannel ( params.genes_dir, "*.all.genes.filtered.fasta", ".all.genes.filtered", { [ [id: it[0].id + '.all.genes', label: 'all', src: it[0].id ], it[1] ] } )
            
        } else {
            // Call genes for all contigs
            ch_gene_call = ch_contigs.map { [ [id: it[0].id + '.all.genes', label: 'all', src: it[0].id ], it[1] ] }
            
            GENE_CALL ( ch_gene_call )
            ch_versions = ch_versions.mix( GENE_CALL.out.versions.first() )

            ch_genes = GENE_CALL.out.genes
        }

        ch_genes_reformatted = ch_genes.map { meta, file -> [ [id: meta.src ], file ]  }

        ch_viral_genes = VIRAL_DETECTION.out.viral_ids.join(ch_genes_reformatted, by:0)
                            .map { [ [id: it[0].id +'.virus.genes', src: it[0].id, label: 'virus'], it[1], it[2] ] }
                            .mix ( 
                                VIRAL_DETECTION.out.plasmid_ids.join(ch_genes_reformatted, by:0)
                                .map { [ [id: it[0].id +'.plasmid.genes', src: it[0].id, label: 'plasmid'], it[1], it[2] ] }
                            )
        
        // ch_viral_genes = ch_genes.join(VIRAL_DETECTION.out.viral_ids, by:0)

        // ch_viral_genes.view()
        // ch_viral_genes = VIRAL_DETECTION.out.viral_ids.map { meta, file -> [ [id: meta.id + '.virus.genes', label: 'virus'], file ]  }
        //                     .mix( VIRAL_DETECTION.out.plasmid_ids.map { meta, file -> [ [id: meta.id + '.plasmid.genes', label: 'plasmid'], file ]  } )

        // ch_viral_genes.view()
        VIRAL_GENE_CALL ( ch_viral_genes )

        VIRAL_GENE_CALL.out.genes.view()
        
        // ch_viral_genes

        // CLUSTER_SEQUENCES ( VIRAL_DETECTION.out.sequences, "mmseqs2" )

        // Build index for abundance estimation (only once)
        // BWA_INDEX ( VIRAL_DETECTION.out.catalogs.concat( ch_genes ) )

        // ch_abundance_input = reads.combine( BWA_INDEX.out.index )
        //                         .map { meta_reads, reads, meta_index, index -> [ [id: meta_reads.id, label: meta_index.id ], reads, index ] }

        // ABUNDANCE ( ch_abundance_input )

        // PROTEIN_CALL ( GENE_CALL.out.gene_catalog )

        // PROTEIN_ANNOTATION ( PROTEIN_CALL.out.protein_catalog, amrfinder_db )

        // VIRAL_ANNOTATION ( VIRAL_DETECTION.out.viral_catalog, virsorter2_db, dram_db, iphop_db )

        // // Collect tables
        // def prepare_table_channels = { preffix, ch ->
        //     ch
        //         .map { [it[1]] }
        //         .collect()
        //         .map { [[id: preffix], it] }
        // }

        // ch_tables = prepare_table_channels('virus', VIRAL_DETECTION.out.virus_tables )
        //             .concat( prepare_table_channels('virus.filtered', VIRAL_DETECTION.out.virus_filtered_tables ) )
        //             .concat( prepare_table_channels('plasmid', VIRAL_DETECTION.out.plasmid_tables ) )
        //             .concat( prepare_table_channels('plasmid.filtered', VIRAL_DETECTION.out.plasmid_filtered_tables ) )
        //             .concat( prepare_table_channels('amg', VIRAL_ANNOTATION.out.amgs ) )
        //             .concat( prepare_table_channels('host.genus', VIRAL_ANNOTATION.out.iphop_genus ) )
        //             .concat( prepare_table_channels('host.genome', VIRAL_ANNOTATION.out.iphop_genomes ) )

        // COLLECT_TABLES ( ch_tables )

        // summary channel versions
        // ch_versions = ASSEMBLY.out.versions
        //                 .mix(VIRAL_DETECTION.out.versions)
        //                 .mix(GENE_CALL.out.versions)
        //                 .mix(ABUNDANCE.out.versions)
        //                 .mix(PROTEIN_CALL.out.versions)
        //                 .mix(PROTEIN_ANNOTATION.out.versions)
        //                 .mix(VIRAL_ANNOTATION.out.versions)


    emit:
        versions = ch_versions

}
