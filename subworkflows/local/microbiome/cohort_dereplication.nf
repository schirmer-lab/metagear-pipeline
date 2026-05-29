include { INPUT_CHECK                  } from "$projectDir/subworkflows/local/common/input_check"
include { ABUNDANCE                    } from "$projectDir/subworkflows/local/common/abundance"

include { STAGE_DREP_WORK              } from "$projectDir/modules/local/cohort/stage_drep_work"
include { DREP_DEREPLICATE             } from "$projectDir/modules/nf-core/drep/dereplicate/main"
include { BUILD_MAG_CATALOG            } from "$projectDir/modules/local/cohort/build_mag_catalog"
include { GATHER_GENOMES_DIR           } from "$projectDir/modules/local/cohort/gather_genomes_dir"
include { GTDBTK_CLASSIFYWF            } from "$projectDir/modules/local/gtdbtk/classifywf/main"
include { SUMMARIZE_MAG_CATALOG        } from "$projectDir/modules/local/cohort/summarize_mag_catalog"

workflow COHORT_DEREPLICATION_INIT {

    main:
        if ( !params.input )    { exit 1, 'Input samplesheet not specified (--input)' }
        if ( !params.bins_dir ) { exit 1, 'Bin directory not specified (--bins_dir)' }
        ch_input = file(params.input)

        gtdb_tk_db = Channel.fromPath("${params.gtdb_tk_db}", checkIfExists: true).first()
        // NOTE: param name `gtdb_tk_db` matches the rest of the pipeline
        // (gene_analysis, MSP). See subworkflows/local/setup/databases.nf.

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


workflow COHORT_DEREPLICATION {

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

        // ─── 2. Cohort dereplication ─────────────────────────────────────────
        // ext.args carries `--S_algorithm skani --genomeInfo drep_work/genomeInfo.csv
        // -comp 50 -con 10` (see conf/metagear/cohort_dereplication.config).
        ch_drep_input = ch_all_bins.map { files -> [ [id: 'cohort'], files ] }
        ch_drep_work  = STAGE_DREP_WORK.out.drep_work.map { dir -> [ [id: 'drep_work'], dir ] }
        DREP_DEREPLICATE ( ch_drep_input, ch_drep_work )
        // DREP emits versions via topic-style channel (not versions.yml).

        // Extract Cdb/Wdb from data_tables/*.csv for downstream joins.
        ch_drep_tables = DREP_DEREPLICATE.out.summary_tables
                                .map { _meta, csvs ->
                                    def csv_list = csvs instanceof List ? csvs : [csvs]
                                    def cdb = csv_list.find { it.name == 'Cdb.csv' }
                                    def wdb = csv_list.find { it.name == 'Wdb.csv' }
                                    return tuple(cdb, wdb)
                                }
        ch_cdb = ch_drep_tables.map { cdb, _wdb -> cdb }
        ch_wdb = ch_drep_tables.map { _cdb, wdb -> wdb }

        // dRep emits dereplicated_genomes/*.fa one path per file. Gather them
        // back into a directory the next steps can consume cleanly.
        ch_drep_reps_list = DREP_DEREPLICATE.out.fastas
                                .map { meta, fas -> [meta, (fas instanceof List ? fas : [fas])] }

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
        // task.ext.extension is set to 'fa' in cohort_dereplication.config.
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
        // Reuse the same genomeInfo CSV that fed dRep (concatenated by
        // STAGE_DREP_WORK).
        ch_genome_info_csv = STAGE_DREP_WORK.out.drep_work
                                .map { dir -> file("${dir}/genomeInfo.csv") }

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
