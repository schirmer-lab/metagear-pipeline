include { INPUT_CHECK                              } from "$projectDir/subworkflows/local/common/input_check"
include { ASSEMBLY                                 } from "$projectDir/subworkflows/local/common/assembly"
include { VIRAL_DETECTION                          } from "$projectDir/subworkflows/local/virus/detection"
include { BACTERIAL_BINNING_INIT; BACTERIAL_BINNING } from "$projectDir/subworkflows/local/microbiome/bacterial_binning"
include { TIARA_TIARA                              } from "$projectDir/modules/nf-core/tiara/tiara/main"
include { EXTRACT_UNBINNED                         } from "$projectDir/modules/local/binette/extract_unbinned/main"
include { MMSEQS_EASYTAXONOMY                      } from "$projectDir/modules/local/mmseqs/easytaxonomy/main"
include { MERGE_CONTIG_CLASSIFICATION              } from "$projectDir/modules/local/metagear/mge/merge_contig_classification"
include { CLASSIFY_GENES                           } from "$projectDir/modules/local/metagear/mge/classify_genes"

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
        // into the dereplication iteration (dRep + GTDB-Tk on the
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
        // Load full-assembly contigs whenever we can:
        //   --contigs_dir set    → load from disk (skip ASSEMBLY)
        //   --chromosome_dir set without contigs_dir → no full assembly available;
        //                          ch_contigs stays empty and Tiara + MERGE fall
        //                          back to ch_chromosome_seqs as their reference.
        //   neither set         → run ASSEMBLY.
        // Having ch_contigs populated even when --chromosome_dir bypasses
        // VIRAL_DETECTION lets Tiara see virus + plasmid contigs (not just the
        // chromosome partition), which is needed for conflict detection in the
        // per-contig TSV (genomad-vs-Tiara disagreements on the viral set).
        ch_contigs = Channel.empty()
        if ( params.contigs_dir ) {
            ch_contigs = createExistingDirChannel ( params.contigs_dir, "*.contigs.fa.gz", ".contigs.fa", false )
        } else if ( !params.chromosome_dir ) {
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
        // unknown}. Run on the full assembly when available so virus + plasmid
        // contigs are also labelled — that's what lets the per-contig TSV
        // surface genomad-vs-Tiara disagreements (e.g. a eukaryotic transposon
        // false-positively flagged as viral by genomad).
        // Falls back to ch_chromosome_seqs only when --chromosome_dir is set
        // without --contigs_dir (no full assembly to use).
        def has_full_contigs = ( !params.chromosome_dir || params.contigs_dir )
        def ch_tiara_input   = has_full_contigs ? ch_contigs : ch_chromosome_seqs
        TIARA_TIARA ( ch_tiara_input )
        ch_tiara_labels = TIARA_TIARA.out.classifications

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

        // ─── 6. MMseqs2 easy-taxonomy on the unbinned contigs (opt-in) ───────
        // Provides contig-level taxonomy for the long-tail of chromosome
        // contigs that didn't reach an MQ+ bin. Off by default — the LCA
        // signal on unbinned contigs is low-confidence (single-protein hits,
        // sub-domain precision discarded by the merge anyway) and the
        // ~200 GiB GTDB DB makes it expensive. Bin-attributable contigs
        // already get trustworthy lineage from GTDB-Tk in dereplication.
        // Enable with --enable_contig_taxonomy true for exploratory analysis.
        def ch_mmseqs_lca = Channel.empty()
        if ( params.enable_contig_taxonomy ) {
            MMSEQS_EASYTAXONOMY ( EXTRACT_UNBINNED.out.unbinned, mmseqs_taxonomy_db )
            ch_versions   = ch_versions.mix( MMSEQS_EASYTAXONOMY.out.versions.first() )
            ch_mmseqs_lca = MMSEQS_EASYTAXONOMY.out.lca
        }

        // ─── 7. Per-contig classification TSV ────────────────────────────────
        // Joins all evidence channels by sample id (String key, since metas
        // drift across subworkflows). Missing channels get assets/empty.txt as
        // the staged fallback; merge_contig_classification.py treats zero-byte
        // / non-directory inputs as empty evidence sets.
        //
        // Reference contigs are the FULL assembly when available (covers
        // viral/plasmid/chromosome all). Only when --chromosome_dir is set
        // WITHOUT --contigs_dir do we lack the full assembly and fall back to
        // the chromosome partition. Same rule as the Tiara routing above —
        // when full contigs are loaded, virus + plasmid rows appear in the
        // TSV and the evidence columns can flag genomad-vs-Tiara conflicts.
        def ch_ref_contigs = has_full_contigs ? ch_contigs : ch_chromosome_seqs

        ch_merge_input = ch_ref_contigs
            .map { meta, fa -> [meta.id, meta, fa] }
            .join( ch_viral_ids.map                          { meta, f -> [meta.id, f] }, by: 0, remainder: true )
            .map { i -> [i[0], i[1], i[2], i[3] ?: empty_file] }
            .join( ch_plasmid_ids.map                        { meta, f -> [meta.id, f] }, by: 0, remainder: true )
            .map { i -> [i[0], i[1], i[2], i[3], i[4] ?: empty_file] }
            .join( BACTERIAL_BINNING.out.bins.map            { meta, files -> [meta.id, files] }, by: 0, remainder: true )
            .map { i -> [i[0], i[1], i[2], i[3], i[4], i[5] ?: empty_file] }
            .join( ch_mmseqs_lca.map                         { meta, f -> [meta.id, f] }, by: 0, remainder: true )
            .map { i -> [i[0], i[1], i[2], i[3], i[4], i[5], i[6] ?: empty_file] }
            .join( ch_tiara_labels.map                       { meta, f -> [meta.id, f] }, by: 0, remainder: true )
            .map { i -> [i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7] ?: empty_file] }
            .map { i -> [i[1], i[2], i[3], i[4], i[5], i[6], i[7]] }   // drop id, keep [meta, contigs, viral, plasmid, bin_dir, mmseqs, tiara]

        MERGE_CONTIG_CLASSIFICATION ( ch_merge_input )
        ch_versions = ch_versions.mix( MERGE_CONTIG_CLASSIFICATION.out.versions.first() )

        // ─── 8. Classify gene-catalog representatives by per-contig class ────
        // Cross-walks the per-contig primary_class onto the gene clusters TSV
        // from gene_analysis. Opt-in: only fires when --gene_clusters_tsv is
        // set (typically auto-injected by the wrapper when a prior gene_analysis
        // run is detected). Produces classification/all.genes.clusters.classified.refined.tsv
        // — a clean re-derivation from per-contig signal; does NOT merge with
        // viral_analysis's `.draft` aggregated TSV. multi_class flags clusters
        // whose members span ≥2 primary_classes (potentially interesting biology
        // like HGT or auxiliary metabolic genes, not necessarily errors).
        ch_classified_genes = Channel.empty()
        if ( params.gene_clusters_tsv ) {
            def gene_clusters_file = file(params.gene_clusters_tsv, checkIfExists: true)
            ch_classify_input = MERGE_CONTIG_CLASSIFICATION.out.tsv
                .map { meta, tsv -> tsv }
                .collect()
                .map { tsvs -> [ [id: 'all.genes'], tsvs, gene_clusters_file ] }

            CLASSIFY_GENES ( ch_classify_input )
            ch_versions         = ch_versions.mix( CLASSIFY_GENES.out.versions )
            ch_classified_genes = CLASSIFY_GENES.out.classified
        }

    emit:
        // Bacterial binning artefacts (per-MAG taxonomy moves to dereplication)
        bins             = BACTERIAL_BINNING.out.bins
        bin_qc_summary   = BACTERIAL_BINNING.out.bin_qc_summary
        unbinned_contigs = EXTRACT_UNBINNED.out.unbinned

        // Eukaryote labels and contig-fallback taxonomy
        tiara_labels         = ch_tiara_labels
        contig_taxonomy_lca  = ch_mmseqs_lca

        // Pass-through for downstream / external consumers
        viral_sequences      = ch_viral_sequences
        plasmid_sequences    = ch_plasmid_sequences
        chromosome_contigs   = ch_chromosome_seqs

        // The integrated_classification headline deliverable: per-contig TSV
        per_contig_tsv       = MERGE_CONTIG_CLASSIFICATION.out.tsv

        // Gene-catalog classification (empty channel unless --gene_clusters_tsv set)
        classified_genes      = ch_classified_genes

        versions = ch_versions
}
