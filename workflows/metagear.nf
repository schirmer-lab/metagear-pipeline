/* --- IMPORT LOCAL SUBWORKFLOWS --- */

include { SETUP } from "$projectDir/workflows/setup"

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { QUALITY_CONTROL_INIT; QUALITY_CONTROL } from "$projectDir/subworkflows/local/common/quality_control"
include { MICROBIAL_PROFILES_INIT; MICROBIAL_PROFILES; METAPHLAN_PROFILES  } from "$projectDir/subworkflows/local/microbiome/microbial_profiles"

include { GENE_ANALYSIS_INIT; GENE_ANALYSIS } from "$projectDir/subworkflows/local/microbiome/gene_analysis"

include { VIRAL_ANALYSIS_INIT; VIRAL_ANALYSIS } from "$projectDir/subworkflows/local/microbiome/viral_analysis"

include { INTEGRATED_CLASSIFICATION_INIT; INTEGRATED_CLASSIFICATION } from "$projectDir/subworkflows/local/microbiome/integrated_classification"

include { COHORT_DEREPLICATION_INIT; COHORT_DEREPLICATION } from "$projectDir/subworkflows/local/microbiome/cohort_dereplication"

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

        if ( params.workflow == "gene_analysis" ) {
            init = GENE_ANALYSIS_INIT ( )

            ch_metaphlan_profiles = Channel.empty()
            if ( init.metaphlan_profiles ) {
                METAPHLAN_PROFILES( init.validated_input, init.metaphlan_db )
                ch_metaphlan_profiles = METAPHLAN_PROFILES.out.merged_profiles.map{ it[1] }
            }

            GENE_ANALYSIS ( init.validated_input, ch_metaphlan_profiles, init.gtdb_tk_db, init.amrfinder_db )
            ch_versions = GENE_ANALYSIS.out.versions
        }

        if ( params.workflow == "viral_analysis" ) {
            init = VIRAL_ANALYSIS_INIT ( )

            VIRAL_ANALYSIS ( init.reads, init.genomad_db, init.checkv_db, init.pharokka_db, init.virsorter2_db, init.dram_db, init.iphop_db, init.amrfinder_db )
            ch_versions = VIRAL_ANALYSIS.out.versions
        }

        // Integrated classification — v1 minimum: assembly + viral/plasmid
        // partitioning (geNomad) + bacterial binning (SemiBin2+MetaBAT2 →
        // Binette → CheckM2 → GTDB-Tk). Contig fallback, plasmid typing,
        // and per-contig TSV merge are deferred to a follow-up.
        if ( params.workflow == "integrated_classification" ) {
            init = INTEGRATED_CLASSIFICATION_INIT ( )

            INTEGRATED_CLASSIFICATION (
                init.reads,
                init.genomad_db,
                init.checkv_db,
                init.checkm2_db,
                init.mmseqs_taxonomy_db,
                init.biome_lookup
            )
            ch_versions = INTEGRATED_CLASSIFICATION.out.versions
        }

        // Cohort dereplication — v2: dRep (skani) on per-sample bins, GTDB-Tk
        // on cluster representatives, coverm-genome MAG×sample abundance, and
        // back-fill of v1's per-contig TSV with the cohort lineage.
        if ( params.workflow == "cohort_dereplication" ) {
            init = COHORT_DEREPLICATION_INIT ( )

            COHORT_DEREPLICATION (
                init.reads,
                init.bins_inputs,
                init.per_contig_tsv,
                init.gtdb_tk_db
            )
            ch_versions = COHORT_DEREPLICATION.out.versions
        }


    emit:
        versions = ch_versions
        summary_data = ch_summary_data

}
