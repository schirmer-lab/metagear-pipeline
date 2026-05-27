include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { createExistingDirChannel; createExistingFileChannel } from "$projectDir/subworkflows/local/utils/existing_data"

include { ASSEMBLY } from "$projectDir/subworkflows/local/common/assembly"

include { VIRAL_DETECTION } from "$projectDir/subworkflows/local/virus/detection"
include { VIRAL_ANNOTATION } from "$projectDir/subworkflows/local/virus/annotation"
include { AMG_POSTPROCESS } from "$projectDir/subworkflows/local/virus/amg_postprocess"

include { GENE_CALL; VIRAL_GENE_CALL } from "$projectDir/subworkflows/local/common/gene_call"
include { PROTEIN_CALL } from "$projectDir/subworkflows/local/common/protein_call"
include { PROTEIN_ANNOTATION } from "$projectDir/subworkflows/local/common/protein_annotation"

include { ABUNDANCE } from "$projectDir/subworkflows/local/common/abundance"

include { CLUSTER_SEQUENCES as CLUSTER_GENES; CLUSTER_SEQUENCES as CLUSTER_PROTEINS; CLUSTER_SEQUENCES as CLUSTER_VIRUS;  CLUSTER_SEQUENCES as CLUSTER_PLASMID } from "$projectDir/subworkflows/local/common/clustering"

include { FIND_REPRESENTATIVES; MERGE_CLUSTER_ANNOTATIONS } from "$projectDir/modules/local/metagear/mge/viral_clusters"

include { COLLECT_TABLES } from "$projectDir/modules/local/metagear/mge/summarize"

workflow VIRAL_ANALYSIS_INIT {

    main:
        if ( params.input ) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

        genomad_db = Channel.fromPath("${params.genomad_db}", checkIfExists: true).first()
        checkv_db = Channel.fromPath("${params.checkv_db}", checkIfExists: true).first()
        pharokka_db = Channel.fromPath("${params.pharokka_db}", checkIfExists: true).first()
        virsorter2_db = Channel.fromPath("${params.virsorter2_db}", checkIfExists: true).first()
        dram_db = Channel.fromPath("${params.dram_db}", checkIfExists: true).first()
        iphop_db = Channel.fromPath("${params.iphop_db}", checkIfExists: true).first()
        amrfinder_db = Channel.fromPath("${params.amrfinder_db}", checkIfExists: true).first()

        INPUT_CHECK ( ch_input, "reads" )

    emit:
        genomad_db
        checkv_db
        pharokka_db
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
        pharokka_db
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
        ch_versions =  ch_versions.mix( VIRAL_DETECTION.out.versions )

        CLUSTER_VIRUS ( VIRAL_DETECTION.out.viral_sequences, "mmseqs2", true )
        ch_versions =  ch_versions.mix( CLUSTER_VIRUS.out.versions )

        CLUSTER_PLASMID ( VIRAL_DETECTION.out.plasmid_sequences, "mmseqs2", true )

        ch_genes = Channel.empty()
        if ( params.genes_dir ) {
            ch_genes = createExistingDirChannel ( params.genes_dir, "*.all.genes.filtered.fasta", ".all.genes.filtered", { [ [id: it[0].id + '.all.genes', label: 'all.genes', src: it[0].id ], it[1] ] } )

        } else {
            // Call genes for all contigs
            ch_gene_call = ch_contigs.map { [ [id: it[0].id + '.all.genes', label: 'all.genes', src: it[0].id ], it[1] ] }

            GENE_CALL ( ch_gene_call )
            ch_versions = ch_versions.mix( GENE_CALL.out.versions.first() )

            ch_genes = GENE_CALL.out.genes
        }

        CLUSTER_GENES ( ch_genes, "mmseqs2", true )
        ch_versions =  ch_versions.mix( CLUSTER_GENES.out.versions )

        PROTEIN_CALL ( CLUSTER_GENES.out.representative )
        ch_versions =  ch_versions.mix(PROTEIN_CALL.out.versions)

        CLUSTER_PROTEINS ( PROTEIN_CALL.out.proteins, "mmseqs2", false )

        PROTEIN_ANNOTATION ( CLUSTER_PROTEINS.out.representative, amrfinder_db )
        ch_versions =  ch_versions.mix(PROTEIN_ANNOTATION.out.versions)

        ch_genes_reformatted = ch_genes.map { meta, file -> [ [id: meta.src ], file ]  }

        // Normalize the viral_ids / plasmid_ids meta to `[id: meta.id]` before
        // joining against ch_genes_reformatted. Nextflow's .join(by: 0) compares
        // the full Map at position 0, and VIRAL_DETECTION.out.{viral,plasmid}_ids
        // inherit extra meta fields from upstream (e.g. `single_end: true` set in
        // input_check.nf for single-end samplesheets, since 446e6c7) that the
        // gene channel doesn't carry. Without this normalization the join
        // produces zero tuples for single-end cohorts and VIRAL_ANNOTATION +
        // AMG_POSTPROCESS silently never fire.
        ch_viral_ids_for_join   = VIRAL_DETECTION.out.viral_ids.map   { meta, ids -> [ [id: meta.id], ids ] }
        ch_plasmid_ids_for_join = VIRAL_DETECTION.out.plasmid_ids.map { meta, ids -> [ [id: meta.id], ids ] }

        ch_viral_genes = ch_viral_ids_for_join.join(ch_genes_reformatted, by:0)
                            .map { [ [id: it[0].id +'.virus.genes', src: it[0].id, label: 'virus'], it[1], it[2] ] }
                            .mix (
                                ch_plasmid_ids_for_join.join(ch_genes_reformatted, by:0)
                                .map { [ [id: it[0].id +'.plasmid.genes', src: it[0].id, label: 'plasmid'], it[1], it[2] ] }
                            )

        VIRAL_GENE_CALL ( ch_viral_genes )
        ch_versions =  ch_versions.mix( VIRAL_GENE_CALL.out.versions )

        ch_viral_representatives = VIRAL_GENE_CALL.out.gene_ids
                .map {meta, file -> [ [id: meta.label + '.genes'], file ]}
                .groupTuple( by:0 )
                .combine (
                    CLUSTER_GENES.out.clusters
                    .join( CLUSTER_GENES.out.representative )
                    .join( PROTEIN_CALL.out.proteins )
                    .map { meta, clusters, representative_genes, representative_proteins -> [ clusters, representative_genes, representative_proteins ] }
                )

        // Get viral/plasmid gene and protein representatives from all clustered genes
        FIND_REPRESENTATIVES ( ch_viral_representatives )

        ch_merge_annotations = FIND_REPRESENTATIVES.out.input_clusters_annotated
                    .map { [it[1]] }
                    .collect()
                    .map { [[id: 'all.genes'], it] }

        MERGE_CLUSTER_ANNOTATIONS ( ch_merge_annotations )

        ch_all_sequences = CLUSTER_GENES.out.representative
                            .concat( CLUSTER_PLASMID.out.representative )
                            .concat( CLUSTER_VIRUS.out.representative )

        ch_abundance_input = reads.combine( ch_all_sequences )
                    .map { meta_reads, reads, meta_sequences, sequences -> [ meta_reads + [label: meta_sequences.id], reads, sequences ] }

        ABUNDANCE ( ch_abundance_input, 'contig', file("$projectDir/assets/empty.txt") )
        ch_versions =  ch_versions.mix( ABUNDANCE.out.versions )

        virus_representative_proteins = FIND_REPRESENTATIVES.out.representative_proteins.filter { meta, _ -> meta.id == 'virus.genes' }.map { it[1] }
        input_viral_annotation = CLUSTER_VIRUS.out.representative.combine( virus_representative_proteins )

        VIRAL_ANNOTATION ( input_viral_annotation, pharokka_db, virsorter2_db, dram_db, iphop_db )
        // ch_versions =  ch_versions.mix(VIRAL_ANNOTATION.out.versions)

        // Collect tables
        def prepare_table_channels = { preffix, ch ->
            ch
                .map { [it[1]] }
                .collect()
                .map { [[id: preffix], it] }
        }

        ch_tables = prepare_table_channels('virus', VIRAL_DETECTION.out.virus_tables )
                    .concat( prepare_table_channels('virus.filtered', VIRAL_DETECTION.out.virus_filtered_tables ) )
                    .concat( prepare_table_channels('plasmid', VIRAL_DETECTION.out.plasmid_tables ) )
                    .concat( prepare_table_channels('plasmid.filtered', VIRAL_DETECTION.out.plasmid_filtered_tables ) )
                    .concat( prepare_table_channels('amg', VIRAL_ANNOTATION.out.amgs ) )
                    .concat( prepare_table_channels('host.genus', VIRAL_ANNOTATION.out.iphop_genus ) )
                    .concat( prepare_table_channels('host.genome', VIRAL_ANNOTATION.out.iphop_genomes ) )

        COLLECT_TABLES ( ch_tables )

        ch_amg_postprocess = COLLECT_TABLES.out.summary.filter { meta, _ -> meta.id == 'amg' }
            .join ( VIRAL_ANNOTATION.out.amg_faa.map { meta, file -> file }.collect().map { [[id: 'amg'], it] } )
            .join ( VIRAL_ANNOTATION.out.amg_fna.map { meta, file -> file }.collect().map { [[id: 'amg'], it] } )
            .join ( ABUNDANCE.out.tpm.filter { meta, _ -> meta.id == 'all.genes_tpm' }.map { [[id: 'amg'], it[1]] } )
            .join ( ABUNDANCE.out.rpkm.filter { meta, _ -> meta.id == 'all.genes_rpkm' }.map { [[id: 'amg'], it[1]] } )
            .join ( ABUNDANCE.out.count.filter { meta, _ -> meta.id == 'all.genes_count' }.map { [[id: 'amg'], it[1]] } )

        // ch_amg_postprocess.view()

        AMG_POSTPROCESS ( ch_amg_postprocess, virus_representative_proteins )


    emit:
        versions = ch_versions

}
