include { INPUT_CHECK                  } from "$projectDir/subworkflows/local/common/input_check"
include { ABUNDANCE                    } from "$projectDir/subworkflows/local/common/abundance"

include { PREPARE_DREP_INPUTS          } from "$projectDir/modules/local/cohort/prepare_drep_inputs"
include { STAGE_DREP_WORK              } from "$projectDir/modules/local/cohort/stage_drep_work"
include { DREP_DEREPLICATE             } from "$projectDir/modules/nf-core/drep/dereplicate/main"
include { BUILD_MAG_CATALOG            } from "$projectDir/modules/local/cohort/build_mag_catalog"
include { GATHER_GENOMES_DIR           } from "$projectDir/modules/local/cohort/gather_genomes_dir"
include { GTDBTK_CLASSIFYWF            } from "$projectDir/modules/local/gtdbtk/classifywf/main"
include { ENRICH_PER_CONTIG_TSV        } from "$projectDir/modules/local/cohort/enrich_per_contig"
include { SUMMARIZE_MAG_CATALOG        } from "$projectDir/modules/local/cohort/summarize_mag_catalog"

include { createExistingDirChannel } from "$projectDir/subworkflows/local/utils/existing_data"


workflow COHORT_DEREPLICATION_INIT {

    main:
        if ( !params.input )    { exit 1, 'Input samplesheet not specified (--input)' }
        if ( !params.bins_dir ) { exit 1, 'Bin directory not specified (--bins_dir)' }
        if ( !params.per_contig_tsv_dir ) {
            exit 1, 'Per-contig TSV directory not specified (--per_contig_tsv_dir)'
        }
        ch_input = file(params.input)

        gtdb_tk_db = Channel.fromPath("${params.gtdb_tk_db}", checkIfExists: true).first()
        // NOTE: param name `gtdb_tk_db` matches the rest of the pipeline
        // (gene_analysis, MSP). See subworkflows/local/setup/databases.nf.

        INPUT_CHECK ( ch_input, "reads" )
        ch_reads = INPUT_CHECK.out.validated_input

        // ─── Per-sample Binette outputs ──────────────────────────────────────
        // bacterial_binning.config publishes BINETTE flat:
        //   ${outdir}/binning/binette/<sample>/final_bins/...
        //   ${outdir}/binning/binette/<sample>/final_bins_quality_reports.tsv
        //   ${outdir}/binning/binette/<sample>/final_contig_to_bin.tsv
        // params.bins_dir points at the `binning/binette/` parent.
        ch_bins_inputs = ch_reads
            .map { meta, _reads ->
                def base = file("${params.bins_dir}/${meta.id}")
                def bins = file("${base}/final_bins")
                def qc   = file("${base}/final_bins_quality_reports.tsv")
                def c2b  = file("${base}/final_contig_to_bin.tsv")
                return [meta, bins, qc, c2b]
            }
            .filter { _meta, bins, qc, _c2b -> bins.isDirectory() && qc.exists() }

        // ─── Per-sample v1 per-contig TSVs ───────────────────────────────────
        // Flat layout from MERGE_CONTIG_CLASSIFICATION's publishDir:
        //   ${outdir}/integrated_classification/per_contig/<sample>.contigs.tsv
        ch_per_contig_tsv = createExistingDirChannel(
            params.per_contig_tsv_dir,
            "*.contigs.tsv",
            ".contigs",
            null
        )

    emit:
        reads          = ch_reads
        bins_inputs    = ch_bins_inputs
        per_contig_tsv = ch_per_contig_tsv
        gtdb_tk_db
        versions       = INPUT_CHECK.out.versions
}


workflow COHORT_DEREPLICATION {

    take:
        reads               // [meta, [r1, r2]]
        bins_inputs         // [meta, bins_dir, qc_tsv, contig_to_bin_tsv]
        per_contig_tsv      // [meta, contigs.tsv]   from v1
        gtdb_tk_db

    main:
        ch_versions = Channel.empty()

        // ─── 1. Per-sample renamer + genomeInfo slice ────────────────────────
        // PREPARE_DREP_INPUTS prefixes each bin with the sample ID and emits
        // the matching genomeInfo CSV slice (genome,completeness,contamination)
        // for dRep's --genomeInfo. Bins dir is staged whole; QC is staged as a
        // file.
        PREPARE_DREP_INPUTS (
            bins_inputs.map { meta, bins_dir, qc, _c2b -> [meta, bins_dir, qc] }
        )

        // Collect ALL renamed bins across samples into a single flat list for
        // dRep (which needs every genome FASTA in one directory).
        ch_all_bins = PREPARE_DREP_INPUTS.out.bins
                            .map { _meta, files -> files }
                            .collect()

        // Concatenate per-sample genomeInfo CSVs into a single cohort-level
        // CSV staged inside drep_work_seed/ for DREP_DEREPLICATE's second
        // input slot.
        ch_per_sample_ginfo = PREPARE_DREP_INPUTS.out.genome_info
                                    .map { _meta, csv -> csv }
                                    .collect()

        STAGE_DREP_WORK ( ch_per_sample_ginfo )

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

        // ─── 5b. Enrich the per-sample contigs.tsv with cohort taxonomy ──────
        // Join per-sample contigs.tsv with per-sample contig_to_bin.tsv first,
        // then layer in the cohort-global Cdb/Wdb + GTDB-Tk summaries.
        ch_enrich_per_sample = per_contig_tsv
                                .map { meta, tsv -> [meta.id, meta, tsv] }
                                .join(
                                    bins_inputs.map { meta, _b, _qc, c2b -> [meta.id, c2b] },
                                    by: 0,
                                    remainder: true
                                )
                                .filter { it[3] != null }
                                .map { _id, meta, tsv, c2b -> [meta, tsv, c2b] }

        ENRICH_PER_CONTIG_TSV (
            ch_enrich_per_sample,
            ch_cdb,
            ch_wdb,
            ch_gtdb_summaries
        )
        ch_versions = ch_versions.mix( ENRICH_PER_CONTIG_TSV.out.versions.first() )

        // ─── 5c. Cohort MAG catalog summary ──────────────────────────────────
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
        per_contig_enriched  = ENRICH_PER_CONTIG_TSV.out.enriched

        // Abundance matrices (MAG × sample)
        abundance_tpm   = ABUNDANCE.out.tpm
        abundance_rpkm  = ABUNDANCE.out.rpkm
        abundance_count = ABUNDANCE.out.count

        versions = ch_versions
}
