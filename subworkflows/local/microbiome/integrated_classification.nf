include { INPUT_CHECK                              } from "$projectDir/subworkflows/local/common/input_check"
include { ASSEMBLY                                 } from "$projectDir/subworkflows/local/common/assembly"
include { VIRAL_DETECTION                          } from "$projectDir/subworkflows/local/virus/detection"
include { BACTERIAL_BINNING_INIT; BACTERIAL_BINNING } from "$projectDir/subworkflows/local/microbiome/bacterial_binning"
include { TIARA_TIARA                              } from "$projectDir/modules/nf-core/tiara/tiara/main"
include { EXTRACT_UNBINNED                         } from "$projectDir/modules/local/binette/extract_unbinned/main"
include { MMSEQS_EASYTAXONOMY                      } from "$projectDir/modules/local/mmseqs/easytaxonomy/main"
include { MERGE_CONTIG_CLASSIFICATION              } from "$projectDir/modules/local/metagear/mge/merge_contig_classification"

include { createExistingDirChannel; createExistingFileChannel } from "$projectDir/subworkflows/local/utils/existing_data"


workflow INTEGRATED_CLASSIFICATION_INIT {

    main:
        if ( !params.input ) { exit 1, 'Input samplesheet not specified!' }
        ch_input = file(params.input)

        // DBs required for the full iter-3 chain:
        //   geNomad + CheckV   → viral/plasmid contig classification
        //   CheckM2            → Binette bin refinement and quality scoring
        //   MMseqs2 taxonomy   → contig-level fallback taxonomy for the
        //                        post-Binette unbinned set (long-tail)
        //
        // GTDB-Tk DB is intentionally NOT loaded — per-MAG taxonomy moves
        // into the cohort_dereplication iteration (dRep + GTDB-Tk on the
        // dereplicated representatives; standard cohort practice).
        genomad_db          = Channel.fromPath("${params.genomad_db}",          checkIfExists: true).first()
        checkv_db           = Channel.fromPath("${params.checkv_db}",           checkIfExists: true).first()
        checkm2_db          = Channel.fromPath("${params.checkm2_db}",          checkIfExists: true).first()
        mmseqs_taxonomy_db  = Channel.fromPath("${params.mmseqs_taxonomy_db}",  checkIfExists: true).first()

        INPUT_CHECK ( ch_input, "reads" )

        // BACTERIAL_BINNING_INIT re-reads the samplesheet for per-sample biome.
        // Kept narrow so meta.biome is NOT injected into INPUT_CHECK's reads
        // channel (handoff §6 trap 2).
        binning_init = BACTERIAL_BINNING_INIT ( )

    emit:
        reads              = INPUT_CHECK.out.validated_input
        genomad_db
        checkv_db
        checkm2_db
        mmseqs_taxonomy_db
        biome_lookup       = binning_init.biome_lookup
        versions           = INPUT_CHECK.out.versions
}


workflow INTEGRATED_CLASSIFICATION {

    take:
        reads               // [meta, [r1, r2]]
        genomad_db
        checkv_db
        checkm2_db
        mmseqs_taxonomy_db
        biome_lookup        // [sample_id, biome]

    main:
        ch_versions = Channel.empty()

        // Fallback file used to satisfy a `path()` input slot when a sample's
        // evidence channel is empty (zero viral IDs, no Tiara labels, etc.).
        // The merge_contig_classification.py script treats zero-byte input
        // files / non-directory bin_dir as empty evidence sets.
        def empty_file = file("$projectDir/assets/empty.txt", checkIfExists: true)

        // ─── 1. Assemble reads → contigs ─────────────────────────────────────
        // --contigs_dir   skips reassembly (uses a prior MEGAHIT run on disk)
        // --chromosome_dir skips both assembly AND VIRAL_DETECTION; in that
        // case ch_contigs is empty and Tiara + MERGE downstream fall back to
        // ch_chromosome_seqs as their reference contigs.
        ch_contigs = Channel.empty()
        if ( params.chromosome_dir ) {
            // ch_contigs stays empty; no consumer downstream of VIRAL_DETECTION needs it.
        } else if ( params.contigs_dir ) {
            ch_contigs = createExistingDirChannel ( params.contigs_dir, "*.contigs.fa.gz", ".contigs.fa", false )
        } else {
            ASSEMBLY ( reads )
            ch_versions = ch_versions.mix( ASSEMBLY.out.versions.first() )
            ch_contigs = ASSEMBLY.out.contigs
        }

        // ─── 2. Viral / plasmid detection → chromosome partition ─────────────
        // VIRAL_DETECTION emits chromosome_sequences = contigs minus the
        // FDR-passed viral + plasmid ID sets (asymmetric policy per handoff §2).
        // When --chromosome_dir bypasses this, viral/plasmid channels are empty
        // unless --viral_ids_dir / --plasmid_ids_dir restore them from disk
        // (typical use: a prior viral_analysis run already published these
        // under results/virus/per_sample/, which the wrapper's auto-reuse
        // picks up).
        ch_chromosome_seqs   = Channel.empty()
        ch_viral_sequences   = Channel.empty()
        ch_plasmid_sequences = Channel.empty()
        ch_viral_ids         = Channel.empty()
        ch_plasmid_ids       = Channel.empty()
        if ( params.chromosome_dir ) {
            ch_chromosome_seqs = createExistingDirChannel (
                params.chromosome_dir,
                "*.chromosome.fna.gz",
                ".chromosome.fna",
                null
            )
        } else {
            VIRAL_DETECTION ( ch_contigs, genomad_db, checkv_db )
            ch_versions          = ch_versions.mix( VIRAL_DETECTION.out.versions )
            ch_chromosome_seqs   = VIRAL_DETECTION.out.chromosome_sequences
            ch_viral_sequences   = VIRAL_DETECTION.out.viral_sequences
            ch_plasmid_sequences = VIRAL_DETECTION.out.plasmid_sequences
            ch_viral_ids         = VIRAL_DETECTION.out.viral_ids
            ch_plasmid_ids       = VIRAL_DETECTION.out.plasmid_ids
        }

        // Explicit override: if the user (or auto-reuse) passed --viral_ids_dir
        // / --plasmid_ids_dir, load those instead. Lets `--chromosome_dir` users
        // still classify virus/plasmid contigs in the per-contig TSV, and lets
        // anyone preserve viral_analysis's MERGE_TABLES outputs across re-runs.
        //
        // Per-sample subdir layout produced by viral_detection.config under
        //   results/virus/per_sample/<sample>/{virus,plasmid}.ids.txt
        // sample id is derived from the parent directory's name, so the
        // helper createExistingDirChannel (which keys off the filename
        // baseName) doesn't fit — using an inline channel here instead.
        if ( params.viral_ids_dir ) {
            ch_viral_ids = Channel.fromPath("${params.viral_ids_dir}/*/virus.ids.txt")
                                .map { f -> [ [id: f.parent.name], f ] }
        }
        if ( params.plasmid_ids_dir ) {
            ch_plasmid_ids = Channel.fromPath("${params.plasmid_ids_dir}/*/plasmid.ids.txt")
                                .map { f -> [ [id: f.parent.name], f ] }
        }

        // ─── 3. Bacterial binning on the chromosome partition ────────────────
        BACTERIAL_BINNING (
            ch_chromosome_seqs,
            reads,
            biome_lookup,
            checkm2_db
        )
        ch_versions = ch_versions.mix( BACTERIAL_BINNING.out.versions )

        // ─── 4. Eukaryote screen (label-only in v1) ──────────────────────────
        // Tiara assigns each contig one of {archaea, bacteria, eukarya, organelle,
        // unknown}. Skipped when --chromosome_dir bypassed assembly+detection
        // (no full-assembly contigs available).
        ch_tiara_labels = Channel.empty()
        if ( !params.chromosome_dir ) {
            TIARA_TIARA ( ch_contigs )
            ch_tiara_labels = TIARA_TIARA.out.classifications
        }

        // ─── 5. Compute the true post-Binette unbinned contigs ───────────────
        // Subtract bin-member contig IDs from the chromosome FASTA. For a
        // sample with no Binette bins, `BACTERIAL_BINNING.out.bins` doesn't
        // emit — we substitute an empty list, and EXTRACT_UNBINNED falls
        // through to "all chromosome contigs are unbinned".
        //
        // Re-key by meta.id (String) because BACTERIAL_BINNING normalizes
        // chromosome's meta against reads' (so a meta-map `.join` would drop
        // samples whose chromosome meta is just `[id:sample]`).
        ch_extract_input = ch_chromosome_seqs
            .map { meta, fa -> [meta.id, meta, fa] }
            .join( BACTERIAL_BINNING.out.bins.map { meta, bins -> [meta.id, bins] }, by: 0, remainder: true )
            .map { items -> [items[1], items[2], items[3] ?: []] }

        EXTRACT_UNBINNED ( ch_extract_input )
        ch_versions = ch_versions.mix( EXTRACT_UNBINNED.out.versions.first() )

        // ─── 6. MMseqs2 easy-taxonomy on the unbinned contigs ────────────────
        // Provides contig-level taxonomy for the long-tail of chromosome
        // contigs that didn't reach an MQ+ bin. EXTRACT_UNBINNED.out.unbinned
        // is optional:true — samples with zero unbinned contigs silently skip
        // this step.
        MMSEQS_EASYTAXONOMY ( EXTRACT_UNBINNED.out.unbinned, mmseqs_taxonomy_db )
        ch_versions = ch_versions.mix( MMSEQS_EASYTAXONOMY.out.versions.first() )

        // ─── 7. Per-contig classification TSV ────────────────────────────────
        // Joins all evidence channels by sample id (String key, since metas
        // drift across subworkflows). Missing channels get assets/empty.txt as
        // the staged fallback; merge_contig_classification.py treats zero-byte
        // / non-directory inputs as empty evidence sets.
        //
        // Reference contigs are the FULL assembly when available (covers
        // viral/plasmid/chromosome all); when --chromosome_dir bypassed
        // assembly+detection, fall back to chromosome contigs only.
        def ch_ref_contigs = params.chromosome_dir ? ch_chromosome_seqs : ch_contigs

        ch_merge_input = ch_ref_contigs
            .map { meta, fa -> [meta.id, meta, fa] }
            .join( ch_viral_ids.map                          { meta, f -> [meta.id, f] }, by: 0, remainder: true )
            .map { i -> [i[0], i[1], i[2], i[3] ?: empty_file] }
            .join( ch_plasmid_ids.map                        { meta, f -> [meta.id, f] }, by: 0, remainder: true )
            .map { i -> [i[0], i[1], i[2], i[3], i[4] ?: empty_file] }
            .join( BACTERIAL_BINNING.out.bins_dir.map        { meta, d -> [meta.id, d] }, by: 0, remainder: true )
            .map { i -> [i[0], i[1], i[2], i[3], i[4], i[5] ?: empty_file] }
            .join( MMSEQS_EASYTAXONOMY.out.lca.map           { meta, f -> [meta.id, f] }, by: 0, remainder: true )
            .map { i -> [i[0], i[1], i[2], i[3], i[4], i[5], i[6] ?: empty_file] }
            .join( ch_tiara_labels.map                       { meta, f -> [meta.id, f] }, by: 0, remainder: true )
            .map { i -> [i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7] ?: empty_file] }
            .map { i -> [i[1], i[2], i[3], i[4], i[5], i[6], i[7]] }   // drop id, keep [meta, contigs, viral, plasmid, bin_dir, mmseqs, tiara]

        MERGE_CONTIG_CLASSIFICATION ( ch_merge_input )
        ch_versions = ch_versions.mix( MERGE_CONTIG_CLASSIFICATION.out.versions.first() )

    emit:
        // Bacterial binning artefacts (per-MAG taxonomy moves to cohort_dereplication)
        bins             = BACTERIAL_BINNING.out.bins
        bins_dir         = BACTERIAL_BINNING.out.bins_dir
        bin_qc_summary   = BACTERIAL_BINNING.out.bin_qc_summary
        unbinned_contigs = EXTRACT_UNBINNED.out.unbinned

        // Eukaryote labels and contig-fallback taxonomy
        tiara_labels         = ch_tiara_labels
        contig_taxonomy_lca  = MMSEQS_EASYTAXONOMY.out.lca

        // Pass-through for downstream / external consumers
        viral_sequences      = ch_viral_sequences
        plasmid_sequences    = ch_plasmid_sequences
        chromosome_contigs   = ch_chromosome_seqs

        // The integrated_classification headline deliverable: per-contig TSV
        per_contig_tsv       = MERGE_CONTIG_CLASSIFICATION.out.tsv

        versions = ch_versions
}
