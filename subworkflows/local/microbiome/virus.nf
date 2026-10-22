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

workflow VIRUS_INIT {

    main:
        if ( params.input ) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

        genomad_db = Channel.fromPath("${params.genomad_db}", checkIfExists: true).first()
        checkv_db = Channel.fromPath("${params.checkv_db}", checkIfExists: true).first()
        pharokka_db = Channel.fromPath("${params.pharokka_db}", checkIfExists: true).first()
        virsorter2_db = Channel.fromPath("${params.virsorter2_db}", checkIfExists: true).first()
        dram_db = Channel.fromPath("${params.dram_db}", checkIfExists: true).first()
        iphop_db = Channel.fromPath("${params.iphop_db}", checkIfExists: true).first()
        amrfinder_db = Channel.fromPath("${params.amrfinder_db}", checkIfExists: true).first()
        phatyp_db = Channel.fromPath("${params.phatyp_db}", checkIfExists: true).first()

        INPUT_CHECK ( ch_input, "reads" )

    emit:
        genomad_db
        checkv_db
        pharokka_db
        virsorter2_db
        dram_db
        iphop_db
        amrfinder_db
        phatyp_db

        reads = INPUT_CHECK.out.validated_input
        versions = INPUT_CHECK.out.versions
}


workflow VIRUS {

    take:
        reads // [meta, reads(fast1, fas2)]
        genomad_db
        checkv_db
        pharokka_db
        virsorter2_db
        dram_db
        iphop_db
        amrfinder_db
        phatyp_db

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

        // Detection is the dominant cost and classification runs it too, so its
        // published output can be reused. Only valid against the same contigs;
        // see docs/workflows/virus.md.
        if ( params.viral_detection_dir ) {
            // Points at a previous run's virus/per_sample; the unfiltered plasmid
            // summary sits beside it under virus/genomad.
            def detection_dir = params.viral_detection_dir
            def detection_root = file ( detection_dir ).parent

            def per_sample_table = { String file_name, Map extra_meta ->
                Channel.fromPath ( "${detection_dir}/*/${file_name}", checkIfExists: true )
                    .map { file -> [ [ id: file.parent.name ] + extra_meta, file ] }
            }

            // label is required by CLUSTER_SEQUENCES.
            ch_viral_sequences         = per_sample_table ( "virus.filtered.fna.gz",        [ label: "virus" ] )
            ch_plasmid_sequences       = per_sample_table ( "plasmid.filtered.fna.gz",      [ label: "plasmid" ] )
            ch_viral_ids               = per_sample_table ( "virus.ids.txt",                [:] )
            ch_plasmid_ids             = per_sample_table ( "plasmid.ids.txt",              [:] )
            ch_virus_tables            = per_sample_table ( "virus.summary.tsv",            [:] )
            ch_virus_filtered_tables   = per_sample_table ( "virus.summary.filtered.tsv",   [:] )
            ch_plasmid_filtered_tables = per_sample_table ( "plasmid.summary.filtered.tsv", [:] )

            // The unfiltered plasmid summary comes straight from GENOMAD_PASS1, which
            // publishes under its own tree instead of per_sample.
            ch_plasmid_tables = Channel.fromPath ( "${detection_root}/genomad/virus/*/*_summary/*_plasmid_summary.tsv", checkIfExists: true )
                .map { file -> [ [ id: file.parent.parent.name ], file ] }

        } else {
            VIRAL_DETECTION ( ch_contigs, genomad_db, checkv_db )
            ch_versions =  ch_versions.mix( VIRAL_DETECTION.out.versions )

            ch_viral_sequences         = VIRAL_DETECTION.out.viral_sequences
            ch_plasmid_sequences       = VIRAL_DETECTION.out.plasmid_sequences
            ch_viral_ids               = VIRAL_DETECTION.out.viral_ids
            ch_plasmid_ids             = VIRAL_DETECTION.out.plasmid_ids
            ch_virus_tables            = VIRAL_DETECTION.out.virus_tables
            ch_virus_filtered_tables   = VIRAL_DETECTION.out.virus_filtered_tables
            ch_plasmid_tables          = VIRAL_DETECTION.out.plasmid_tables
            ch_plasmid_filtered_tables = VIRAL_DETECTION.out.plasmid_filtered_tables
        }

        // Default is 'vclust', which implements the MIUViG species criterion the
        // catalog intends: 95% ANI over 85% of the shorter sequence, on merged
        // local alignments. Set --virus_clustering_method mmseqs2 to reproduce
        // catalogs built before this changed; see docs/workflows/virus.md.
        CLUSTER_VIRUS ( ch_viral_sequences, params.virus_clustering_method ?: "vclust", true )
        ch_versions =  ch_versions.mix( CLUSTER_VIRUS.out.versions )

        CLUSTER_PLASMID ( ch_plasmid_sequences, "mmseqs2", true )

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

        ch_representative_genes = Channel.empty()
        ch_gene_clusters = Channel.empty()

        // Re-clustering picks different representatives and silently invalidates the
        // classification and MSP tables. Needs the cluster map too, not just the reps.
        if ( params.representative_genes && params.gene_clusters_tsv ) {
            ch_representative_genes = createExistingFileChannel ( params.representative_genes, { [ [id: "all.genes"], it ] } )
            ch_gene_clusters = createExistingFileChannel ( params.gene_clusters_tsv, { [ [id: "all.genes"], it ] } )

        } else {
            CLUSTER_GENES ( ch_genes, "mmseqs2", true )
            ch_versions =  ch_versions.mix( CLUSTER_GENES.out.versions )

            ch_representative_genes = CLUSTER_GENES.out.representative
            ch_gene_clusters = CLUSTER_GENES.out.clusters
        }

        // Always runs: FIND_REPRESENTATIVES needs the unclustered translation, and
        // that is not a reusable artifact.
        PROTEIN_CALL ( ch_representative_genes )
        ch_versions =  ch_versions.mix(PROTEIN_CALL.out.versions)

        ch_representative_proteins = Channel.empty()

        if ( params.representative_proteins ) {
            ch_representative_proteins = createExistingFileChannel ( params.representative_proteins, { [ [id: "all.genes"], it ] } )

        } else {
            CLUSTER_PROTEINS ( PROTEIN_CALL.out.proteins, "mmseqs2", false )
            ch_representative_proteins = CLUSTER_PROTEINS.out.representative
        }

        if ( !params.representative_proteins_annotations ) {
            PROTEIN_ANNOTATION ( ch_representative_proteins, amrfinder_db )
            ch_versions =  ch_versions.mix(PROTEIN_ANNOTATION.out.versions)
        }

        ch_genes_reformatted = ch_genes.map { meta, file -> [ [id: meta.src ], file ]  }

        // Normalize the viral_ids / plasmid_ids meta to `[id: meta.id]` before
        // joining against ch_genes_reformatted. Nextflow's .join(by: 0) compares
        // the full Map at position 0, and VIRAL_DETECTION.out.{viral,plasmid}_ids
        // inherit extra meta fields from upstream (e.g. `single_end: true` set in
        // input_check.nf for single-end samplesheets, since 446e6c7) that the
        // gene channel doesn't carry. Without this normalization the join
        // produces zero tuples for single-end cohorts and VIRAL_ANNOTATION +
        // AMG_POSTPROCESS silently never fire.
        ch_viral_ids_for_join   = ch_viral_ids.map   { meta, ids -> [ [id: meta.id], ids ] }
        ch_plasmid_ids_for_join = ch_plasmid_ids.map { meta, ids -> [ [id: meta.id], ids ] }

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
                    ch_gene_clusters
                    .join( ch_representative_genes )
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

        ch_all_sequences = ch_representative_genes
                            .concat( CLUSTER_PLASMID.out.representative )
                            .concat( CLUSTER_VIRUS.out.representative )

        ch_abundance_input = reads.combine( ch_all_sequences )
                    .map { meta_reads, reads, meta_sequences, sequences -> [ meta_reads + [label: meta_sequences.id], reads, sequences ] }

        ABUNDANCE ( ch_abundance_input, 'contig', file("$projectDir/assets/empty.txt"), null )
        ch_versions =  ch_versions.mix( ABUNDANCE.out.versions )

        virus_representative_proteins = FIND_REPRESENTATIVES.out.representative_proteins.filter { meta, _ -> meta.id == 'virus.genes' }.map { it[1] }
        input_viral_annotation = CLUSTER_VIRUS.out.representative.combine( virus_representative_proteins )

        VIRAL_ANNOTATION ( input_viral_annotation, pharokka_db, virsorter2_db, dram_db, iphop_db, phatyp_db )
        // ch_versions =  ch_versions.mix(VIRAL_ANNOTATION.out.versions)

        // Collect tables
        def prepare_table_channels = { preffix, ch ->
            ch
                .map { [it[1]] }
                .collect()
                .map { [[id: preffix], it] }
        }

        ch_tables = prepare_table_channels('virus', ch_virus_tables )
                    .concat( prepare_table_channels('virus.filtered', ch_virus_filtered_tables ) )
                    .concat( prepare_table_channels('plasmid', ch_plasmid_tables ) )
                    .concat( prepare_table_channels('plasmid.filtered', ch_plasmid_filtered_tables ) )
                    .concat( prepare_table_channels('amg', VIRAL_ANNOTATION.out.amgs ) )
                    .concat( prepare_table_channels('host.genus', VIRAL_ANNOTATION.out.iphop_genus ) )
                    .concat( prepare_table_channels('host.genome', VIRAL_ANNOTATION.out.iphop_genomes ) )
                    // PhaTYP lifestyle calls take the same route as the host
                    // predictions: per-chunk tables merged by COLLECT_TABLES into
                    // one catalog-level table, so no bespoke merge step is needed.
                    .concat( prepare_table_channels('lifestyle', VIRAL_ANNOTATION.out.lifestyle ) )

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
