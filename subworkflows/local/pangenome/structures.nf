include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"
include { createExistingFileChannel } from "$projectDir/subworkflows/local/utils/existing_data"

include { SEQKIT_SPLIT2 as SEQKIT_SPLIT2_STRUCTURES } from "$projectDir/modules/nf-core/seqkit/split2"
include { SEQKIT_SPLIT2 as SEQKIT_SPLIT2_COMPARE    } from "$projectDir/modules/nf-core/seqkit/split2"

include { FILTER_STRUCTURES_INPUTS } from "$projectDir/modules/local/metagear/utils/filter_structures_inputs"
include { PACKAGE_PHOLD_OFFSITE    } from "$projectDir/modules/local/metagear/utils/package_phold_offsite"
include { PHOLD_PREDICT            } from "$projectDir/modules/local/phold/predict/main"
include { MERGE_PHOLD_PREDICTIONS  } from "$projectDir/modules/local/metagear/utils/merge_phold_predictions"
include { SPLIT_3DI_BY_AA          } from "$projectDir/modules/local/metagear/utils/split_3di_by_aa"
include { PHOLD_COMPARE            } from "$projectDir/modules/local/phold/compare/main"
include { MERGE_PHOLD_COMPARE      } from "$projectDir/modules/local/metagear/utils/merge_phold_compare"
include { MERGE_VIRAL_PHOLD        } from "$projectDir/modules/local/metagear/utils/merge_viral_phold"


workflow STRUCTURES_INIT {

    main:
        if ( !params.input )                              { exit 1, 'Input samplesheet not specified (--input)' }
        if ( !params.representative_proteins )            { exit 1, 'Representative-proteins catalog not specified (--representative_proteins)' }
        if ( !params.representative_proteins_clusters )   { exit 1, 'Protein cluster TSV not specified (--representative_proteins_clusters)' }
        if ( !params.viral_representative_proteins )      { exit 1, 'Viral representative-proteins catalog not specified (--viral_representative_proteins)' }

        // phold_db is only strictly required when this run actually invokes
        // PHOLD (i.e. NOT a --structures_prepare_for_gpu run, where the DB
        // lives on the remote GPU machine). For the resume path
        // (--phold_predict_dir set, no PHOLD_PREDICT on the source side),
        // PHOLD_COMPARE still needs it. So we require phold_db unless
        // structures_prepare_for_gpu is true.
        if ( !params.phold_db && !params.structures_prepare_for_gpu ) {
            exit 1, 'PHOLD structural DB not specified (--phold_db). Install via `metagear download_databases --databases structures` and point phold_db at the resulting directory in ~/.metagear/metagear.config.'
        }

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

        // phold_db channel — present unless we're in --structures_prepare_for_gpu
        // mode (where the DB lives on the remote machine instead).
        phold_db = params.phold_db \
            ? Channel.fromPath("${params.phold_db}", checkIfExists: true).first() \
            : Channel.empty()

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
        phold_db          // path             — PHOLD structural DB directory (or empty when prepare_for_gpu)

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

        // Initialise downstream channels as empty; populated below depending
        // on which branch we take.
        ch_offsite_bundle = Channel.empty()
        ch_all_phold      = Channel.empty()
        ch_viral_phold    = Channel.empty()
        ch_di_fasta       = Channel.empty()

        // ─── Branch A: prepare offsite GPU bundle and stop ───────────────────
        // When the user can't run the GPU step locally (e.g. cluster without
        // Singularity, no SLURM-from-job), this branch packages everything
        // PHOLD_PREDICT needs into one rsync-friendly directory and stops.
        // The user runs the bundled script on a GPU machine, rsyncs the
        // results back, then re-invokes `metagear structures` without
        // --structures_prepare_for_gpu — auto-reuse picks up the rsynced
        // outputs and the pipeline resumes from MERGE_PHOLD_PREDICTIONS
        // onward (Branch B).
        if ( params.structures_prepare_for_gpu ) {

            // Collect all shards into one tuple for PACKAGE_PHOLD_OFFSITE.
            ch_shards_collected = SEQKIT_SPLIT2_STRUCTURES.out.reads
                                    .map { _meta, gz ->
                                        (gz instanceof java.nio.file.Path) ? [gz] : (gz as List)
                                    }
                                    .flatten()
                                    .collect()

            // Join shards with the viral_join_table so PACKAGE has both.
            // FILTER_STRUCTURES_INPUTS.out.viral_join_table is per-meta;
            // we re-key to a sortable channel here.
            ch_package_in = ch_shards_collected
                                .map { shards -> tuple( [id: "cohort"], shards ) }
                                .join( FILTER_STRUCTURES_INPUTS.out.viral_join_table )

            // Templates live in assets/structures/offsite/. README + runner +
            // sbatch wrapper are rendered with cohort-specific values (shard
            // count, ETA estimates per GPU class, suggested walltime).
            readme_template = file("$projectDir/assets/structures/offsite/README.md.in")
            script_template = file("$projectDir/assets/structures/offsite/run_phold_predict.sh.in")
            sbatch_template = file("$projectDir/assets/structures/offsite/submit_phold_predict.sbatch.in")

            PACKAGE_PHOLD_OFFSITE ( ch_package_in, readme_template, script_template, sbatch_template )
            ch_versions       = ch_versions.mix( PACKAGE_PHOLD_OFFSITE.out.versions )
            ch_offsite_bundle = PACKAGE_PHOLD_OFFSITE.out.bundle

        } else {

            // ─── Branch B: PHOLD_PREDICT — either offsite-rsynced or run-here ─
            ch_merge_predict_in = Channel.empty()

            if ( params.phold_predict_dir ) {

                // Resume path: per-shard predict_<id>/ dirs are on disk
                // already (rsynced back from the GPU server). Skip
                // PHOLD_PREDICT entirely and stream the existing dirs
                // straight into MERGE_PHOLD_PREDICTIONS. The merge module
                // tolerates missing files per shard with warnings, so a
                // partial rsync surfaces as a warning rather than a silent
                // skip — and the downstream PHOLD_COMPARE will fail loudly
                // if the merged predict_dir is incomplete.
                ch_merge_predict_in = Channel.fromPath(
                                            "${params.phold_predict_dir}/predict_*",
                                            type: 'dir',
                                            checkIfExists: true
                                        )
                                        .collect()
                                        .map { dirs -> tuple( [id: "cohort"], dirs ) }

            } else {

                // Normal path: run PHOLD_PREDICT per shard (scatter).
                // phold_db is passed as a value channel — PHOLD validates
                // the DB on startup (even for predict, which only uses
                // ProstT5 weights from the same dir), so every shard needs
                // it bound.
                ch_predict_in = SEQKIT_SPLIT2_STRUCTURES.out.reads
                                    .flatMap { meta, gz ->
                                        def files = (gz instanceof java.nio.file.Path) ? [gz] : (gz as List)
                                        files.collect { f ->
                                            def fn = f.getFileName().toString()
                                            def chunkId = fn.replaceFirst(/\.faa\.gz$/, '').replaceFirst(/\.fa\.gz$/, '').replaceFirst(/\.fasta\.gz$/, '')
                                            tuple( [id: chunkId, src: meta.id], f )
                                        }
                                    }

                PHOLD_PREDICT ( ch_predict_in, phold_db )
                ch_versions = ch_versions.mix( PHOLD_PREDICT.out.versions.first() )

                ch_merge_predict_in = PHOLD_PREDICT.out.predict_dir
                                        .map { _meta, dir -> dir }
                                        .collect()
                                        .map { dirs -> tuple( [id: "cohort"], dirs ) }
            }

            // ─── Merge per-shard prediction dirs → one merged predict_dir ────
            MERGE_PHOLD_PREDICTIONS ( ch_merge_predict_in )
            ch_versions = ch_versions.mix( MERGE_PHOLD_PREDICTIONS.out.versions )

            // ─── PHOLD_COMPARE scatter-gather ────────────────────────────────
            // Foldseek search at sensitivity 9.5 scales linearly with input
            // protein count, but each `phold proteins-compare` invocation
            // pays a fixed ~30-60 s DB-load overhead. Sweet spot is
            // 50-200 k proteins per chunk; default 100 k via
            // `structures_compare_shard_size`.
            //
            // For the 911 k-protein mouse cohort this is ~10 chunks; for a
            // 500-sample, ~50 M-protein cohort it's ~500 chunks running
            // 30-60 min each, parallelised across a fat node.

            //  Step a: split cohort AA FASTA into M chunks (chunk size set via
            //          ext.args in conf/metagear/structures.config). The
            //          nf-core SEQKIT_SPLIT2 module branches on
            //          `meta.single_end`; tag it true so we hit the
            //          single-file path instead of the paired-end path that
            //          tries to pass `--read2 null`. Matches the trick the
            //          PREDICT scatter (SEQKIT_SPLIT2_STRUCTURES) already uses.
            ch_split_compare_in = FILTER_STRUCTURES_INPUTS.out.subset_fasta
                                    .map { meta, fa -> tuple( meta + [single_end: true], fa ) }
            SEQKIT_SPLIT2_COMPARE ( ch_split_compare_in )
            ch_versions = ch_versions.mix( SEQKIT_SPLIT2_COMPARE.out.versions )

            //  Step b: per chunk, build a chunk-scoped predict_dir from the
            //          cohort merged_predict_dir by filtering to the chunk's
            //          protein IDs.
            ch_aa_chunks = SEQKIT_SPLIT2_COMPARE.out.reads
                                .flatMap { meta, chunks ->
                                    def cs = (chunks instanceof java.nio.file.Path) ? [chunks] : (chunks as List)
                                    cs.collect { c ->
                                        def fn = c.getFileName().toString()
                                        def chunkId = fn.replaceFirst(/\.faa\.gz$/, '')
                                                       .replaceFirst(/\.fa\.gz$/, '')
                                                       .replaceFirst(/\.fasta\.gz$/, '')
                                                       .replaceFirst(/\.faa$/, '')
                                                       .replaceFirst(/\.fa$/, '')
                                        tuple( [id: chunkId, src: meta.id], c )
                                    }
                                }
            ch_split_3di_in = ch_aa_chunks
                                .combine( MERGE_PHOLD_PREDICTIONS.out.merged_predict_dir
                                            .map { _meta, d -> d } )
            SPLIT_3DI_BY_AA ( ch_split_3di_in )
            ch_versions = ch_versions.mix( SPLIT_3DI_BY_AA.out.versions.first() )

            //  Step c: PHOLD_COMPARE per chunk (Foldseek search).
            ch_compare_in = ch_aa_chunks
                                .join( SPLIT_3DI_BY_AA.out.chunk_predict_dir )
            PHOLD_COMPARE ( ch_compare_in, phold_db )
            ch_versions = ch_versions.mix( PHOLD_COMPARE.out.versions.first() )

            //  Step d: gather per-chunk per_cds_tsv into the cohort TSV.
            ch_merge_compare_in = PHOLD_COMPARE.out.per_cds_tsv
                                    .map { _meta, tsv -> tsv }
                                    .collect()
                                    .map { tsvs -> tuple( [id: "cohort"], tsvs ) }
            MERGE_PHOLD_COMPARE ( ch_merge_compare_in )
            ch_versions  = ch_versions.mix( MERGE_PHOLD_COMPARE.out.versions )
            ch_all_phold = MERGE_PHOLD_COMPARE.out.per_cds_tsv

            //  di_fasta: take from MERGE_PHOLD_PREDICTIONS — its phold_3di.fasta
            //  already covers all cohort proteins (no need to gather per-chunk
            //  copies, which would be duplicates of the merged file's contents).
            ch_di_fasta = MERGE_PHOLD_PREDICTIONS.out.di_fasta

            // ─── Build virus.proteins.phold.tsv via the viral join table ─────
            ch_merge_viral_in = ch_all_phold
                                .join( FILTER_STRUCTURES_INPUTS.out.viral_join_table )
            MERGE_VIRAL_PHOLD ( ch_merge_viral_in )
            ch_versions    = ch_versions.mix( MERGE_VIRAL_PHOLD.out.versions )
            ch_viral_phold = MERGE_VIRAL_PHOLD.out.viral_phold
        }

    emit:
        offsite_bundle  = ch_offsite_bundle        // [meta, offsite_predict/] (only when prepare_for_gpu)
        all_phold       = ch_all_phold             // [meta, all.proteins.phold.tsv] (else)
        viral_phold     = ch_viral_phold           // [meta, virus.proteins.phold.tsv] (else)
        di_fasta        = ch_di_fasta              // 3Di FASTA — reusable downstream (else)
        versions        = ch_versions
}
