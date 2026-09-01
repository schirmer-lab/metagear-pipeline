include { BWA_INDEX                          } from "$projectDir/modules/nf-core/bwa/index/main"
include { COVERM_MAKE                        } from "$projectDir/modules/local/coverm/make"
include { SAMTOOLS_INDEX                     } from "$projectDir/modules/nf-core/samtools/index/main"
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS } from "$projectDir/modules/nf-core/metabat2/jgisummarizebamcontigdepths/main"
include { METABAT2_METABAT2                  } from "$projectDir/modules/nf-core/metabat2/metabat2/main"
include { SEMIBIN_SINGLEEASYBIN              } from "$projectDir/modules/nf-core/semibin/singleeasybin/main"
include { BINETTE                            } from "$projectDir/modules/local/binette/main"


workflow BACTERIAL_BINNING_INIT {

    main:
        // Re-read here rather than adding `biome` to INPUT_CHECK's meta, which would
        // invalidate every meta-keyed cache downstream of input_check.nf.
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
        // Reuse channels carry only `id` while reads' meta carries `single_end` too,
        // so joining on the maps drops every sample. Re-key both by meta.id.
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
        // SAMTOOLS_INDEX, METABAT2_METABAT2 and SEMIBIN_SINGLEEASYBIN are
        // topic-versions only; reading .out.versions on them aborts the DAG.
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
        // biome joins the meta only here, so no other meta-keyed task sees the change.
        // --environment is wired via the ext.args2 closure in the config.
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
        // Binette takes --bin_dirs XOR --contig2bin_tables, never both.
        // SemiBin's meta carries `biome` and the others' does not, so re-key by
        // meta.id again. remainder:true keeps a sample whose binner produced 0 bins.
        ch_binette_input = ch_chromosome
            .map { meta, fa -> [ meta.id, meta, fa ] }
            .join( SEMIBIN_SINGLEEASYBIN.out.output_fasta.map { meta, bins -> [ meta.id, bins ] }, by: 0, remainder: true )
            .join( METABAT2_METABAT2.out.fasta.map            { meta, bins -> [ meta.id, bins ] }, by: 0, remainder: true )
            .map { items -> [ items[1], items[2], items[3] ?: [], items[4] ?: [] ] }  // [meta, contigs, semibin_bins_or_empty, metabat_bins_or_empty]

        BINETTE ( ch_binette_input, checkm2_db )
        ch_versions = ch_versions.mix( BINETTE.out.versions.first() )

        // NOTE: GTDB-Tk taxonomy is intentionally NOT run here. Per-sample
        // GTDB-Tk on per-sample bins is wasteful at cohort scale (10–30× more
        // runs than needed) and produces redundant calls because dRep at 95%
        // ANI is by construction the GTDB species threshold — every bin in a
        // dRep cluster maps to the same species. The canonical taxonomy step
        // moves into the mag iteration: dRep dereplicates
        // across samples, then GTDB-Tk runs once on the representatives.

    emit:
        bins             = BINETTE.out.bins
        bin_qc_summary   = BINETTE.out.quality
        // The TRUE post-Binette unbinned set is computed in
        // classification.nf via EXTRACT_UNBINNED (chromosome FASTA
        // minus all bin-member contig IDs). MetaBAT2's `.unbinned` is just
        // its own loose-bin set, not what the per-contig TSV needs.
        versions         = ch_versions
}
