include { BWA_INDEX                          } from "$projectDir/modules/nf-core/bwa/index/main"
include { COVERM_MAKE                        } from "$projectDir/modules/local/coverm/make"
include { SAMTOOLS_INDEX                     } from "$projectDir/modules/nf-core/samtools/index/main"
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS } from "$projectDir/modules/nf-core/metabat2/jgisummarizebamcontigdepths/main"
include { METABAT2_METABAT2                  } from "$projectDir/modules/nf-core/metabat2/metabat2/main"
include { SEMIBIN_SINGLEEASYBIN              } from "$projectDir/modules/nf-core/semibin/singleeasybin/main"
include { BINETTE                            } from "$projectDir/modules/local/binette/main"


workflow BACTERIAL_BINNING_INIT {

    main:
        // Per integrated_classification_handoff §6 trap 2: do not pull `biome` out
        // of INPUT_CHECK's meta — that would invalidate every meta-keyed cache
        // downstream of input_check.nf. Instead, re-read the samplesheet here and
        // join the biome onto our own channels before SemiBin2.
        if ( !params.input ) { exit 1, 'Input samplesheet not specified!' }

        ch_biome_lookup = Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [ row.sample, (row.biome?.trim() ?: 'global') ] }

    emit:
        biome_lookup = ch_biome_lookup
}


workflow BACTERIAL_BINNING {

    take:
        chromosome_contigs   // [meta, fasta]            meta.id = sample id
        reads                // [meta, [r1, r2]]
        biome_lookup         // [sample_id, biome]       from _INIT
        checkm2_db

    main:
        ch_versions = Channel.empty()

        // ─── 0. Normalize chromosome_contigs' meta against reads' meta ─────
        // When chromosome_contigs comes from createExistingDirChannel (i.e.
        // --contigs_dir or --chromosome_dir reuse), its meta is reconstructed
        // from the filename and only contains `id`. Reads' meta carries
        // `single_end` (and whatever INPUT_CHECK injected). The two are NOT
        // equal as Groovy maps, so `.join(by:0)` would silently drop every
        // sample. Re-key both by meta.id (String), pick reads' meta, then
        // use ch_chromosome for the rest of the subworkflow.
        ch_chromosome = reads
            .map { meta, _r -> [ meta.id, meta ] }
            .join( chromosome_contigs.map { meta, fa -> [ meta.id, fa ] }, by: 0 )
            .map { _id, meta, fa -> [ meta, fa ] }

        // ─── 1. Index sample's chromosome contigs (per-sample) ──────────────
        BWA_INDEX ( ch_chromosome )
        ch_versions = ch_versions.mix( BWA_INDEX.out.versions.first() )

        // ─── 2. Map sample's own reads to its own contigs (single-sample) ───
        //         No cross-sample mapping — see handoff §4 (decided).
        ch_reads_with_index = reads
            .join( BWA_INDEX.out.index, by: 0 )       // [meta, reads, bwa_idx_dir]

        COVERM_MAKE ( ch_reads_with_index, true )      // true = ref_is_index
        ch_versions = ch_versions.mix( COVERM_MAKE.out.versions.first() )

        // ─── 3. Index the BAM (BAI required by metabat2 depth step) ─────────
        // SAMTOOLS_INDEX uses topic-style version emit (no versions.yml path) —
        // not mixed into ch_versions here. Same for METABAT2_METABAT2 and
        // SEMIBIN_SINGLEEASYBIN below. Topic-version collection is a separate
        // pipeline-level concern; not blocking this subworkflow.
        SAMTOOLS_INDEX ( COVERM_MAKE.out.alignments )

        // ─── 4. Depth summary per contig for MetaBAT2 ───────────────────────
        ch_bam_bai = COVERM_MAKE.out.alignments
            .join( SAMTOOLS_INDEX.out.index, by: 0 )  // [meta, bam, bai/csi/crai]

        METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( ch_bam_bai )
        ch_versions = ch_versions.mix( METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions.first() )

        // ─── 5. MetaBAT2 binning ────────────────────────────────────────────
        ch_metabat_input = ch_chromosome
            .join( METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth, by: 0 )  // [meta, fasta, depth]

        METABAT2_METABAT2 ( ch_metabat_input )
        // METABAT2_METABAT2 uses topic-style versions — see SAMTOOLS_INDEX note above.

        // ─── 6. SemiBin2 binning (biome-aware via meta.biome) ───────────────
        // Inject biome into meta JUST BEFORE SemiBin2 so the meta change is
        // scoped to this process (no other meta-keyed task downstream sees it).
        // SEMIBIN_SINGLEEASYBIN's tag uses meta.id only; --environment is wired
        // via ext.args2 closure in conf/metagear/bacterial_binning.config.
        ch_semibin_input = ch_chromosome
            .join( COVERM_MAKE.out.alignments, by: 0 )         // [meta, fasta, bam]
            .map { meta, fasta, bam -> [ meta.id, meta, fasta, bam ] }
            .join( biome_lookup, by: 0 )                        // [id, meta, fasta, bam, biome]
            .map { id, meta, fasta, bam, biome ->
                [ meta + [ biome: biome ], fasta, bam ]
            }

        SEMIBIN_SINGLEEASYBIN ( ch_semibin_input )
        // SEMIBIN_SINGLEEASYBIN uses topic-style versions — see SAMTOOLS_INDEX note above.

        // ─── 7. Binette refinement ─────────────────────────────────────────
        //   Binette 1.2.1's CLI requires --bin_dirs XOR --contig2bin_tables
        //   (not both — verified at runtime: "Either ... must be provided, but
        //   not both"). Both upstream binners emit per-bin FASTAs in
        //   directories, so we use --bin_dirs for both and stage SemiBin's
        //   output_bins/*.fa.gz next to MetaBAT2's bin FASTAs.
        //
        //   --min_completeness 50 (set via ext.args in config) gates to MIMAG
        //   MQ+ bins, so final_bins/ is already the passing set.
        //
        //   Meta normalization: SEMIBIN_SINGLEEASYBIN's input meta carries
        //   `biome` (so its ext.args2 closure can read meta.biome to set
        //   --environment), and that enriched meta propagates to its outputs.
        //   ch_chromosome and METABAT2_METABAT2.out.fasta still have the slim
        //   [id, single_end] meta. Re-key all three by meta.id (String) so the
        //   join is meta-shape-agnostic; ch_chromosome's clean meta becomes
        //   the canonical meta for BINETTE.
        ch_binette_input = ch_chromosome
            .map { meta, fa -> [ meta.id, meta, fa ] }
            .join( SEMIBIN_SINGLEEASYBIN.out.output_fasta.map { meta, bins -> [ meta.id, bins ] }, by: 0 )
            .join( METABAT2_METABAT2.out.fasta.map            { meta, bins -> [ meta.id, bins ] }, by: 0 )
            .map { _id, meta, fa, semibin_bins, metabat_bins -> [ meta, fa, semibin_bins, metabat_bins ] }  // [meta, contigs, semibin_bins, metabat_bins]

        BINETTE ( ch_binette_input, checkm2_db )
        ch_versions = ch_versions.mix( BINETTE.out.versions.first() )

        // NOTE: GTDB-Tk taxonomy is intentionally NOT run here. Per-sample
        // GTDB-Tk on per-sample bins is wasteful at cohort scale (10–30× more
        // runs than needed) and produces redundant calls because dRep at 95%
        // ANI is by construction the GTDB species threshold — every bin in a
        // dRep cluster maps to the same species. The canonical taxonomy step
        // moves into the cohort_dereplication iteration: dRep dereplicates
        // across samples, then GTDB-Tk runs once on the representatives.

    emit:
        bins             = BINETTE.out.bins
        bin_qc_summary   = BINETTE.out.quality
        // The TRUE post-Binette unbinned set is computed in
        // integrated_classification.nf via EXTRACT_UNBINNED (chromosome FASTA
        // minus all bin-member contig IDs). MetaBAT2's `.unbinned` is just
        // its own loose-bin set, not what the per-contig TSV needs.
        versions         = ch_versions
}
