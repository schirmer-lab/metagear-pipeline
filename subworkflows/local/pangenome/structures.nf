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

        // Not required only when preparing a GPU bundle: the DB lives on the remote machine.
        if ( !params.phold_db && !params.structures_prepare_for_gpu ) {
            exit 1, 'PHOLD structural DB not specified (--phold_db). Install via `metagear download_databases --databases structures` and point phold_db at the resulting directory in ~/.metagear/metagear.config.'
        }

        ch_input = file(params.input)

        all_proteins   = createExistingFileChannel ( params.representative_proteins,         { [ [id: "cohort"], it ] } )
        clusters_tsv   = createExistingFileChannel ( params.representative_proteins_clusters, { [ [id: "cohort"], it ] } )
        viral_proteins = createExistingFileChannel ( params.viral_representative_proteins,    { [ [id: "cohort"], it ] } )

        // Optional: without it every rep counts as unannotated, so those scopes collapse to `all`.
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

        // Scope selection and viral top-up both live in bin/filter_structures_inputs.py.
        ch_filter_in = all_proteins
                            .join( viral_proteins )
                            .join( clusters_tsv )
                            .join( pfam_tsv )
        FILTER_STRUCTURES_INPUTS ( ch_filter_in )
        ch_versions = ch_versions.mix( FILTER_STRUCTURES_INPUTS.out.versions )

        // Shard for parallel ProstT5 inference; chunk size, not file count, sets the shard count.
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

        // Package what PHOLD_PREDICT needs and stop; re-invoking without the flag resumes.
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

                // Shards already on disk from the GPU server: skip PHOLD_PREDICT.
                ch_merge_predict_in = Channel.fromPath(
                                            "${params.phold_predict_dir}/predict_*",
                                            type: 'dir',
                                            checkIfExists: true
                                        )
                                        .collect()
                                        .map { dirs -> tuple( [id: "cohort"], dirs ) }

            } else {

                // phold_db is a value channel: PHOLD validates it on startup, even for predict.
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

            // Each compare pays a fixed DB-load cost, so chunk at ~100k proteins rather than per file.

            // single_end:true picks SEQKIT_SPLIT2's single-file path; the paired one passes --read2 null.
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
