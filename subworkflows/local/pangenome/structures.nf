include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"
include { createExistingFileChannel } from "$projectDir/subworkflows/local/utils/existing_data"

include { SEQKIT_SPLIT2 as SEQKIT_SPLIT2_STRUCTURES } from "$projectDir/modules/nf-core/seqkit/split2"

include { FILTER_STRUCTURES_INPUTS } from "$projectDir/modules/local/metagear/utils/filter_structures_inputs"
include { PHOLD_PREDICT            } from "$projectDir/modules/local/phold/predict/main"
include { MERGE_PHOLD_PREDICTIONS  } from "$projectDir/modules/local/metagear/utils/merge_phold_predictions"
include { PHOLD_COMPARE            } from "$projectDir/modules/local/phold/compare/main"
include { MERGE_VIRAL_PHOLD        } from "$projectDir/modules/local/metagear/utils/merge_viral_phold"


workflow STRUCTURES_INIT {

    main:
        if ( !params.input )                              { exit 1, 'Input samplesheet not specified (--input)' }
        if ( !params.representative_proteins )            { exit 1, 'Representative-proteins catalog not specified (--representative_proteins)' }
        if ( !params.representative_proteins_clusters )   { exit 1, 'Protein cluster TSV not specified (--representative_proteins_clusters)' }
        if ( !params.viral_representative_proteins )      { exit 1, 'Viral representative-proteins catalog not specified (--viral_representative_proteins)' }
        if ( !params.phold_db )                      { exit 1, 'PHOLD structural DB not specified (--phold_db). Install via `metagear download_databases --databases structures` and point phold_db at the resulting directory in ~/.metagear/metagear.config.' }

        ch_input = file(params.input)

        all_proteins   = createExistingFileChannel ( params.representative_proteins,         { [ [id: "cohort"], it ] } )
        clusters_tsv   = createExistingFileChannel ( params.representative_proteins_clusters, { [ [id: "cohort"], it ] } )
        viral_proteins = createExistingFileChannel ( params.viral_representative_proteins,    { [ [id: "cohort"], it ] } )

        // Pfam annotations are optional — the filter degrades gracefully:
        // when missing, every all-protein-rep is treated as unannotated, so
        // `unannotated` and `unannotated_plus_duf` collapse to `all`. We
        // surface this as an empty channel so the join in STRUCTURES still
        // composes; the FILTER step handles the null on disk.
        pfam_tsv = params.representative_proteins_annotations \
            ? createExistingFileChannel ( params.representative_proteins_annotations, { [ [id: "cohort"], it ] } ) \
            : Channel.of ( [ [id: "cohort"], file("$projectDir/assets/empty.txt") ] )

        phold_db = Channel.fromPath("${params.phold_db}", checkIfExists: true).first()

        INPUT_CHECK ( ch_input, "reads" )

    emit:
        all_proteins
        clusters_tsv
        viral_proteins
        pfam_tsv
        phold_db
        versions = INPUT_CHECK.out.versions
}


workflow STRUCTURES {

    take:
        all_proteins      // [meta, fa.gz]    — all.proteins.representative.fa.gz
        clusters_tsv      // [meta, tsv]      — all.proteins.clusters.tsv
        viral_proteins    // [meta, fa.gz]    — virus.proteins.representative.fa.gz
        pfam_tsv          // [meta, tsv]      — all.proteins.FG_IPS_Pfam.tsv (or empty)
        phold_db     // path             — PHOLD structural DB directory

    main:
        ch_versions = Channel.empty()

        // ─── 1. Pick the subset to feed PHOLD + build viral join table ───────
        // The scope selector (params.structures_scope) plus the viral top-up
        // logic both live in bin/filter_structures_inputs.py. Output is one
        // FASTA destined for PHOLD_PREDICT and a mapping TSV used at the end
        // for the viral catalog table.
        ch_filter_in = all_proteins
                            .join( viral_proteins )
                            .join( clusters_tsv )
                            .join( pfam_tsv )
        FILTER_STRUCTURES_INPUTS ( ch_filter_in )
        ch_versions = ch_versions.mix( FILTER_STRUCTURES_INPUTS.out.versions )

        // ─── 2. Shard the subset for parallel ProstT5 prediction ─────────────
        // ProstT5 inference is the heavy step. We split into N chunks via
        // SEQKIT_SPLIT2 (configurable through structures.config's ext.args
        // -s <chunk-size>) so PHOLD_PREDICT runs N-way parallel on multi-core
        // / multi-GPU hosts. The shard count is gated by chunk size, not
        // file count — same pattern as iPHoP.
        ch_split_in = FILTER_STRUCTURES_INPUTS.out.subset_fasta
                            .map { meta, fa -> tuple( meta + [single_end: true], fa ) }
        SEQKIT_SPLIT2_STRUCTURES ( ch_split_in )
        ch_versions = ch_versions.mix( SEQKIT_SPLIT2_STRUCTURES.out.versions )

        // Per-shard inputs to PHOLD_PREDICT. SEQKIT_SPLIT2 emits a (meta, [fa, fa, ...])
        // tuple; flatMap unfolds into one tuple per shard, tagging each shard
        // with its own meta.id (taken from the chunk filename minus `.fa.gz`).
        ch_predict_in = SEQKIT_SPLIT2_STRUCTURES.out.reads
                            .flatMap { meta, gz ->
                                def files = (gz instanceof java.nio.file.Path) ? [gz] : (gz as List)
                                files.collect { f ->
                                    def fn = f.getFileName().toString()
                                    def chunkId = fn.replaceFirst(/\.faa\.gz$/, '').replaceFirst(/\.fa\.gz$/, '').replaceFirst(/\.fasta\.gz$/, '')
                                    tuple( [id: chunkId, src: meta.id], f )
                                }
                            }

        // ─── 3. PHOLD_PREDICT per shard (scatter) ────────────────────────────
        // phold_db is passed as a value channel — PHOLD validates the DB on
        // startup (even for predict, which only uses ProstT5 weights from
        // the same dir), so every shard needs it bound.
        PHOLD_PREDICT ( ch_predict_in, phold_db )
        ch_versions = ch_versions.mix( PHOLD_PREDICT.out.versions.first() )

        // ─── 4. Collect per-shard prediction dirs → one merged predict_dir ───
        ch_merge_predict_in = PHOLD_PREDICT.out.predict_dir
                            .map { _meta, dir -> dir }
                            .collect()
                            .map { dirs -> tuple( [id: "cohort"], dirs ) }
        MERGE_PHOLD_PREDICTIONS ( ch_merge_predict_in )
        ch_versions = ch_versions.mix( MERGE_PHOLD_PREDICTIONS.out.versions )

        // ─── 5. PHOLD_COMPARE (Foldseek vs PHOLD DB) — single cohort job ─────
        // Joins the merged 3Di predictions with the original AA FASTA from
        // FILTER_STRUCTURES_INPUTS so PHOLD can compute coverage / identity
        // metrics against the original sequences.
        ch_compare_in = FILTER_STRUCTURES_INPUTS.out.subset_fasta
                            .join( MERGE_PHOLD_PREDICTIONS.out.merged_predict_dir )
        PHOLD_COMPARE ( ch_compare_in, phold_db )
        ch_versions = ch_versions.mix( PHOLD_COMPARE.out.versions )

        // ─── 6. Build virus.proteins.phold.tsv via the viral join table ──────
        ch_merge_viral_in = PHOLD_COMPARE.out.per_cds_tsv
                            .join( FILTER_STRUCTURES_INPUTS.out.viral_join_table )
        MERGE_VIRAL_PHOLD ( ch_merge_viral_in )
        ch_versions = ch_versions.mix( MERGE_VIRAL_PHOLD.out.versions )

    emit:
        all_phold       = PHOLD_COMPARE.out.per_cds_tsv     // [meta, all.proteins.phold.tsv]
        viral_phold     = MERGE_VIRAL_PHOLD.out.viral_phold // [meta, virus.proteins.phold.tsv]
        di_fasta        = PHOLD_COMPARE.out.di_fasta        // 3Di FASTA — reusable downstream
        versions        = ch_versions
}
