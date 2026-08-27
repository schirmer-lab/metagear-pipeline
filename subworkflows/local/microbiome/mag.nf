include { INPUT_CHECK                  } from "$projectDir/subworkflows/local/common/input_check"
include { ABUNDANCE                    } from "$projectDir/subworkflows/local/common/abundance"

include { STAGE_DREP_WORK              } from "$projectDir/modules/local/cohort/stage_drep_work"
include { SKANI_TRIANGLE               } from "$projectDir/modules/local/cohort/skani_triangle"
include { BATCH_BINS                   } from "$projectDir/modules/local/cohort/batch_bins"
include { DREP_DEREPLICATE             } from "$projectDir/modules/nf-core/drep/dereplicate/main"
include { MERGE_DREP_TABLES            } from "$projectDir/modules/local/cohort/merge_drep_tables"
include { BUILD_MAG_CATALOG            } from "$projectDir/modules/local/cohort/build_mag_catalog"
include { GATHER_GENOMES_DIR           } from "$projectDir/modules/local/cohort/gather_genomes_dir"
include { GTDBTK_CLASSIFYWF            } from "$projectDir/modules/local/gtdbtk/classifywf/main"
include { SUMMARIZE_MAG_CATALOG        } from "$projectDir/modules/local/cohort/summarize_mag_catalog"

workflow MAG_INIT {

    main:
        if ( !params.input )    { exit 1, 'Input samplesheet not specified (--input)' }
        if ( !params.bins_dir ) { exit 1, 'Bin directory not specified (--bins_dir)' }
        ch_input = file(params.input)

        gtdb_tk_db = Channel.fromPath("${params.gtdb_tk_db}", checkIfExists: true).first()
        // NOTE: param name `gtdb_tk_db` matches the rest of the pipeline
        // (genes, MSP). See subworkflows/local/setup/databases.nf.

        INPUT_CHECK ( ch_input, "reads" )
        ch_reads = INPUT_CHECK.out.validated_input

        // ─── Per-sample Binette outputs ──────────────────────────────────────
        // bacterial_binning.config publishes BINETTE's per-sample deliverables
        // (renamed and flattened) to:
        //   ${outdir}/assemblies/bins/<sample>/
        //     ├── <sample>_binN.fa                  (MAGs, --prefix-ed by Binette)
        //     ├── <sample>.quality_report.tsv       (renamed at publish time;
        //     │                                      sample-prefixed so cohort
        //     │                                      collection doesn't collide)
        //     └── <sample>.contig_to_bin.tsv        (renamed at publish time;
        //                                            sample-prefixed for the
        //                                            same reason)
        // params.bins_dir points at the `assemblies/bins/` parent. The
        // sample subdir holds the .fa files at top level — no final_bins/
        // nesting — so the cohort step below globs them directly.
        // The c2b path (sample contig_to_bin.tsv) was previously surfaced
        // here for ENRICH_PER_CONTIG_TSV; that step was dropped (see comment
        // in the main workflow below) so we no longer need it in the tuple.
        ch_bins_inputs = ch_reads
            .map { meta, _reads ->
                def sample_root = file("${params.bins_dir}/${meta.id}")
                def qc          = file("${sample_root}/${meta.id}.quality_report.tsv")
                return [meta, sample_root, qc]
            }
            .filter { _meta, bins, qc -> bins.isDirectory() && qc.exists() }

    emit:
        reads          = ch_reads
        bins_inputs    = ch_bins_inputs
        gtdb_tk_db
        versions       = INPUT_CHECK.out.versions
}


workflow MAG {

    take:
        reads               // [meta, [r1, r2]]
        bins_inputs         // [meta, bins_dir, qc_tsv]
        gtdb_tk_db

    main:
        ch_versions = Channel.empty()

        // ─── 1. Collect cohort bins + cohort genomeInfo CSV ──────────────────
        // Bins are globally unique already (BINETTE was invoked with
        // `--prefix <sample>`), so no per-sample rename step is needed —
        // PREPARE_DREP_INPUTS used to handle that, plus a per-sample QC TSV →
        // genomeInfo CSV slice. Both have been absorbed: bins are globbed
        // directly from the published sample dirs, and STAGE_DREP_WORK now
        // does the CSV transformation in one cohort-level pass.
        ch_all_bins = bins_inputs
                            .map { _meta, sample_root, _qc -> files("${sample_root}/*.fa") }
                            .flatten()
                            .collect()

        ch_all_qc = bins_inputs
                            .map { _meta, _root, qc -> qc }
                            .collect()

        STAGE_DREP_WORK ( ch_all_qc )

        // ─── 2a. Cohort pairwise ANI (skani triangle) ────────────────────────
        // skani triangle bulk-compares every bin against every other bin and
        // emits a sparse TSV of pairs above its internal ANI floor. This is
        // O(N²) but ~1000× cheaper per pair than dRep's full pipeline, so it
        // scales to tens of thousands of bins where a single all-vs-all dRep
        // would OOM or time out (see resources.yaml comment on DREP_DEREPLICATE).
        ch_skani_in = ch_all_bins.map { files -> [ [id: 'cohort'], files ] }
        SKANI_TRIANGLE ( ch_skani_in )
        ch_versions = ch_versions.mix( SKANI_TRIANGLE.out.versions )

        // ─── 2b. Partition bins into dRep batches ────────────────────────────
        // Connected components at params.drep_batch_ani (default 0.90)
        // strictly contain dRep's 0.95 secondary clustering threshold, so
        // within-batch dRep produces the same cohort-wide representatives
        // as a single all-vs-all run on the full cohort. See bin/batch_bins.py
        // for the full safety argument. Singleton bins become 1-bin batches
        // — they still run through dRep so its -comp/-con quality filter is
        // applied uniformly across the cohort.
        ch_genome_info_csv = STAGE_DREP_WORK.out.drep_work
                                .map { dir -> file("${dir}/genomeInfo.csv") }

        // Two-input call: the first channel zips (ani, genome_info) by meta;
        // the second is the cohort bin list (a single-emission list channel).
        // We deliberately do NOT use .combine() to thread the bins through:
        // .combine() unpacks list-typed emissions into cartesian factors, so
        // a list of 326 paths becomes 326 separate emissions — which would
        // run BATCH_BINS 326 times with 1 bin each. Passing as a separate
        // input keeps the list intact as a single emission.
        ch_batch_bins_in = SKANI_TRIANGLE.out.ani
                                .join( ch_genome_info_csv.map { csv -> [ [id: 'cohort'], csv ] } )
        BATCH_BINS ( ch_batch_bins_in, ch_all_bins )
        ch_versions = ch_versions.mix( BATCH_BINS.out.versions )

        // ─── 2c. Per-batch dRep (scatter) ────────────────────────────────────
        // Flatten the glob of batch_<NNN>/ dirs into one emission per batch,
        // then split into the two channels DREP_DEREPLICATE expects (bins,
        // drep_work seed). The two channels emit in lockstep so dRep's
        // process pairing matches them by position.
        ch_per_batch = BATCH_BINS.out.batches
                            .map { meta, batch_dirs ->
                                (batch_dirs instanceof List ? batch_dirs : [batch_dirs])
                            }
                            .flatten()
                            .map { dir ->
                                def bid  = dir.name           // batch_<NNN>
                                def bins = file("${dir}/bins/*.fa")
                                def seed = file("${dir}/drep_work_seed")
                                tuple( [id: bid], bins, seed )
                            }

        ch_drep_input = ch_per_batch.map { meta, bins, _seed -> tuple(meta, bins) }
        ch_drep_work  = ch_per_batch.map { meta, _bins, seed -> tuple([id: "${meta.id}_work"], seed) }

        DREP_DEREPLICATE ( ch_drep_input, ch_drep_work )
        ch_versions = ch_versions.mix(DREP_DEREPLICATE.out.versions)
        // DREP emits versions via topic-style channel (not versions.yml).

        // ─── 2d. Concatenate per-batch Cdb/Wdb into cohort tables (gather) ───
        // dRep emits batch-local cluster IDs in Cdb.csv / Wdb.csv. We rename
        // each batch's CSVs to `batch_<NNN>.Cdb.csv` / `batch_<NNN>.Wdb.csv`
        // via collectFile so MERGE_DREP_TABLES can re-prefix the per-batch
        // rows with a `batch_id` column (keeps cluster IDs disambiguable in
        // the cohort table).
        ch_drep_tables_renamed = DREP_DEREPLICATE.out.summary_tables
                                .map { meta, csvs ->
                                    def csv_list = csvs instanceof List ? csvs : [csvs]
                                    def cdb = csv_list.find { it.name == 'Cdb.csv' }
                                    def wdb = csv_list.find { it.name == 'Wdb.csv' }
                                    return tuple(meta.id, cdb, wdb)
                                }

        ch_cdb_renamed = ch_drep_tables_renamed
                                .collectFile() { bid, cdb, _wdb -> [ "${bid}.Cdb.csv", cdb.text ] }
        ch_wdb_renamed = ch_drep_tables_renamed
                                .collectFile() { bid, _cdb, wdb -> [ "${bid}.Wdb.csv", wdb.text ] }

        // Two-input call (same pattern as BATCH_BINS): .combine() spreads
        // list-typed emissions into cartesian factors, so combining two
        // .collect()'d lists doesn't produce a single (cdbs, wdbs) tuple
        // — it flattens. Passing as separate inputs lets Nextflow zip them
        // by position and keeps each list intact as one emission.
        ch_cdb_collected = ch_cdb_renamed.collect().map { cdbs -> tuple([id: 'cohort'], cdbs) }
        MERGE_DREP_TABLES ( ch_cdb_collected, ch_wdb_renamed.collect() )
        ch_versions = ch_versions.mix( MERGE_DREP_TABLES.out.versions )

        ch_cdb = MERGE_DREP_TABLES.out.cdb.map { _meta, cdb -> cdb }
        ch_wdb = MERGE_DREP_TABLES.out.wdb.map { _meta, wdb -> wdb }

        // ─── 2e. Gather per-batch winners into a single cohort channel ───────
        // dRep emits dereplicated_genomes/*.fa one path per file per batch.
        // Flatten across batches and re-emit as a single meta-tagged tuple
        // that the rest of the workflow can consume unchanged.
        ch_drep_reps_list = DREP_DEREPLICATE.out.fastas
                                .map { _meta, fas -> (fas instanceof List ? fas : [fas]) }
                                .flatten()
                                .collect()
                                .map { fas -> tuple([id: 'cohort'], fas) }

        // ─── 3. Build the cohort MAG catalog + contig→MAG lookup ─────────────
        // BUILD_MAG_CATALOG accepts a flat list of FASTAs and emits
        //   mag_catalog.fa.gz, contig_to_mag.tsv, mag_contig_lengths.tsv
        ch_reps_flat = ch_drep_reps_list.map { _meta, fas -> fas }.flatten().collect()
        BUILD_MAG_CATALOG ( ch_reps_flat )
        ch_versions = ch_versions.mix( BUILD_MAG_CATALOG.out.versions )

        // ─── 4. Cohort abundance (genome mode) ───────────────────────────────
        // ABUNDANCE in 'genome' mode runs coverm genome with
        // --genome-definition contig_to_mag.tsv, producing per-MAG matrices.
        // BWA_INDEX + COVERM_MAKE are shared with contig mode, so existing
        // callers' caches remain intact (see commit log).
        ch_abundance_in = reads.combine( BUILD_MAG_CATALOG.out.catalog )
                            .map { meta_reads, reads_fastqs, catalog ->
                                // Fixed cohort label; ABUNDANCE uses it as
                                // the BWA index key so it's built once.
                                [ meta_reads + [label: 'cohort_mag_catalog'], reads_fastqs, catalog ]
                            }

        ABUNDANCE (
            ch_abundance_in,
            'genome',
            BUILD_MAG_CATALOG.out.contig_to_mag
        )
        ch_versions = ch_versions.mix( ABUNDANCE.out.versions )

        // ─── 5a. GTDB-Tk on cluster representatives ──────────────────────────
        // dRep emits winners as <sample>.binette_binN.fa; classifywf's
        // task.ext.extension is set to 'fa' in mag.config.
        GATHER_GENOMES_DIR ( ch_drep_reps_list )

        ch_gtdbtk_in = GATHER_GENOMES_DIR.out.dir.combine( gtdb_tk_db )
                            .map { meta, dir, db -> [meta, dir, db] }
        GTDBTK_CLASSIFYWF ( ch_gtdbtk_in, false )
        ch_versions = ch_versions.mix( GTDBTK_CLASSIFYWF.out.versions.first() )

        // Collect the per-domain summary TSVs (bac120, ar53). GTDB-Tk emits one
        // or both depending on what was in the cohort.
        ch_gtdb_summaries = GTDBTK_CLASSIFYWF.out.summary
                                .map { _meta, summary ->
                                    summary instanceof List ? summary : [summary]
                                }
                                .flatten()
                                .collect()

        // Note: a per-sample per-contig TSV enrichment step (joining cohort
        // Cdb/Wdb + GTDB-Tk lineages back into each sample's contigs.tsv)
        // used to live here. It was dropped — the join is a few-line pandas
        // operation that analysts can do on demand from the published
        // building blocks (Cdb, Wdb, GTDB summary, per-sample contig_to_bin,
        // per_contig TSV). See docs/recipes/per-contig-taxonomy-enrichment.md.

        // ─── 5b. Cohort MAG catalog summary ──────────────────────────────────
        // Reuses ch_genome_info_csv from STAGE_DREP_WORK (declared above in
        // section 2b — same channel can feed multiple consumers in DSL2).
        SUMMARIZE_MAG_CATALOG (
            ch_cdb,
            ch_wdb,
            ch_genome_info_csv,
            ch_gtdb_summaries
        )
        ch_versions = ch_versions.mix( SUMMARIZE_MAG_CATALOG.out.versions )

    emit:
        // Cohort headline deliverables
        mag_catalog          = BUILD_MAG_CATALOG.out.catalog
        mag_catalog_summary  = SUMMARIZE_MAG_CATALOG.out.catalog_summary
        cluster_table        = ch_cdb
        cluster_winners      = ch_wdb
        gtdbtk_summary       = ch_gtdb_summaries

        // Abundance matrices (MAG × sample)
        abundance_tpm   = ABUNDANCE.out.tpm
        abundance_rpkm  = ABUNDANCE.out.rpkm
        abundance_count = ABUNDANCE.out.count

        versions = ch_versions
}
