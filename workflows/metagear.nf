/* --- IMPORT LOCAL SUBWORKFLOWS --- */

include { SETUP } from "$projectDir/workflows/setup"

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { QUALITY_CONTROL_INIT; QUALITY_CONTROL } from "$projectDir/subworkflows/local/common/quality_control"
include { MICROBIAL_PROFILES_INIT; MICROBIAL_PROFILES; METAPHLAN_PROFILES  } from "$projectDir/subworkflows/local/microbiome/microbial_profiles"

include { GENES_INIT; GENES } from "$projectDir/subworkflows/local/microbiome/genes"

include { VIRUS_INIT; VIRUS } from "$projectDir/subworkflows/local/microbiome/virus"

include { CLASSIFICATION_INIT; CLASSIFICATION } from "$projectDir/subworkflows/local/microbiome/classification"

include { MAG_INIT; MAG } from "$projectDir/subworkflows/local/microbiome/mag"

include { MSP_INIT; MSP } from "$projectDir/subworkflows/local/pangenome/msp"

include { STRUCTURES_INIT; STRUCTURES } from "$projectDir/subworkflows/local/pangenome/structures"

/* --- RUN MAIN WORKFLOW --- */
workflow METAGEAR {

    main:

        ch_versions = Channel.empty()
        ch_summary_data = Channel.empty()

        if (params.workflow == null || params.workflow.trim().isEmpty()) {
            INPUT_CHECK ( file(params.input), "reads" )
        }

        // Setup handler
        if ( params.workflow == "download_databases" ) {
            SETUP ( )
            ch_versions = SETUP.out.versions
        }

        // Quality Control handler
        if ( params.workflow.startsWith("qc_") ) {
            init = QUALITY_CONTROL_INIT ( )
            QUALITY_CONTROL ( init.validated_input, init.kneaddata_refdb )
            ch_versions = QUALITY_CONTROL.out.versions

            ch_summary_data = QUALITY_CONTROL.out.fastqc_zip_pre.collect{it[1]}.ifEmpty([])
                    .mix(QUALITY_CONTROL.out.fastqc_zip_post.collect{it[1]}.ifEmpty([]))
                    .mix(QUALITY_CONTROL.out.summary_plot.collect{it}.ifEmpty([]))
        }

        // Microbial profiles
        if ( params.workflow == "microbial_profiles" ) {
            init = MICROBIAL_PROFILES_INIT ( )
            MICROBIAL_PROFILES ( init.validated_input, init.metaphlan_db, init.uniref90_db, init.chocoplhan_db )
            ch_versions = MICROBIAL_PROFILES.out.versions
        }

        if ( params.workflow == "genes" ) {
            init = GENES_INIT ( )

            GENES ( init.validated_input, init.amrfinder_db )
            ch_versions = GENES.out.versions
        }

        // MSP — gene-centric species (MetaSpecies Pangenomes). MSPminer
        // co-abundance clustering on the gene catalog, GTDB-Tk taxonomy on
        // MSP representative sequences, and MetaPhlAn cross-walk. Builds on a
        // prior `genes` run; reads the gene catalog and per-sample abundance
        // matrices from disk.
        if ( params.workflow == "msp" ) {
            init = MSP_INIT ( )

            ch_metaphlan_profiles = Channel.empty()
            if ( init.metaphlan_profiles ) {
                METAPHLAN_PROFILES( init.validated_input, init.metaphlan_db )
                ch_metaphlan_profiles = METAPHLAN_PROFILES.out.merged_profiles.map{ it[1] }
            }

            MSP (
                init.representative_genes,
                init.representative_genes_count,
                init.representative_genes_rpkm,
                init.gtdb_tk_db,
                ch_metaphlan_profiles
            )
            ch_versions = MSP.out.versions
        }

        if ( params.workflow == "virus" ) {
            init = VIRUS_INIT ( )

            VIRUS ( init.reads, init.genomad_db, init.checkv_db, init.pharokka_db, init.virsorter2_db, init.dram_db, init.iphop_db, init.amrfinder_db )
            ch_versions = VIRUS.out.versions
        }

        // Integrated classification — v1 minimum: assembly + viral/plasmid
        // partitioning (geNomad) + bacterial binning (SemiBin2+MetaBAT2 →
        // Binette → CheckM2 → GTDB-Tk). Contig fallback, plasmid typing,
        // and per-contig TSV merge are deferred to a follow-up.
        if ( params.workflow == "classification" ) {
            init = CLASSIFICATION_INIT ( )

            CLASSIFICATION (
                init.reads,
                init.genomad_db,
                init.checkv_db,
                init.checkm2_db,
                init.mmseqs_taxonomy_db,
                init.biome_lookup
            )
            ch_versions = CLASSIFICATION.out.versions
        }

        // Dereplication — v2: dRep (skani) on per-sample bins, GTDB-Tk
        // on cluster representatives, coverm-genome MAG×sample abundance, and
        // a cohort MAG catalog summary. Per-contig taxonomy enrichment is left
        // as a downstream pandas join — see docs/recipes/per-contig-taxonomy-enrichment.md.
        if ( params.workflow == "mag" ) {
            init = MAG_INIT ( )

            MAG (
                init.reads,
                init.bins_inputs,
                init.gtdb_tk_db
            )
            ch_versions = MAG.out.versions
        }

        // Structures — protein structural-homology annotation via PHOLD
        // (ProstT5 → Foldseek). Builds on a prior virus/genes run; reads
        // the cohort protein catalog + Pfam table + protein clusters +
        // viral catalog from disk. Emits per-rep PHOLD annotations for
        // the full catalog and the viral catalog separately.
        if ( params.workflow == "structures" ) {
            init = STRUCTURES_INIT ( )

            STRUCTURES (
                init.all_proteins,
                init.clusters_tsv,
                init.viral_proteins,
                init.pfam_tsv,
                init.phold_db
            )
            ch_versions = STRUCTURES.out.versions
        }


    emit:
        versions = ch_versions
        summary_data = ch_summary_data

}
