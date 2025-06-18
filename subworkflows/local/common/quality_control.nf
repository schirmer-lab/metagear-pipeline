/* --- Perform quality control on reads --- */

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"
include { CONCATENATE_READS } from "$projectDir/modules/local/metagear/utils/concatenate_reads"

include { FASTQC as FASTQC_RAW; FASTQC as FASTQC_CLEAN} from "$projectDir/modules/nf-core/fastqc/main"
include { TRIMGALORE } from "$projectDir/modules/nf-core/trimgalore/main"
include { KNEADDATA; PARSE_KNEADDATA; SUMMARY_KNEADDATA } from "$projectDir/modules/local/kneaddata/main"
include { FIX_FASTX_HEADERS } from "$projectDir/modules/local/metagear/utils/fix_fastx_headers"

/* --- Main Workflow --- */
workflow QUALITY_CONTROL_INIT {
    main:

        if ( params.input ) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

        INPUT_CHECK ( ch_input, "reads" )

        // Get the kneaddata reference database
        def kneaddata_dbs = params.kneaddata_refdb.collect { file(it) }
        kneaddata_refdb = Channel.value(kneaddata_dbs)

        CONCATENATE_READS ( INPUT_CHECK.out.validated_input )

    emit:
        validated_input = CONCATENATE_READS.out
        kneaddata_refdb = kneaddata_refdb
        versions = INPUT_CHECK.out.versions

}

/* --- Main Workflow --- */
workflow QUALITY_CONTROL {

    take:
        validated_input // meta, reads
        kneaddata_refdb

    main:
        // FastQC, quality trimming, adapter removal, and host decontamination

        FASTQC_RAW ( validated_input )

        TRIMGALORE ( validated_input )

        // Fix fastq headers (enabled by default)
        if ( params.fix_fastq_header ) {
            FIX_HEADER ( TRIMGALORE.out.reads )
            ch_kneaddata = FIX_HEADER.out.reads
        }else {
            ch_kneaddata = TRIMGALORE.out.reads
        }

        KNEADDATA ( ch_kneaddata , kneaddata_refdb )

        PARSE_KNEADDATA ( KNEADDATA.out.kneaddata_log )

        PARSE_KNEADDATA.out.kneadata_stats
            .map { meta, kneaddata_stats -> kneaddata_stats }
            .collect()
            .set { all_kneaddata_stats }

        SUMMARY_KNEADDATA ( all_kneaddata_stats )

        // Append "_clean" to IDs for post cleaning fastqc
        ch_fastqc_clean = KNEADDATA.out.reads
                            .map { meta, reads -> [[id: meta.id + "_clean"] , reads] }

        FASTQC_CLEAN ( ch_fastqc_clean )

        ch_versions_qc = FASTQC_RAW.out.versions
                            .mix(TRIMGALORE.out.versions)
                            .mix(KNEADDATA.out.versions)
                            .mix(FASTQC_CLEAN.out.versions)

    emit:
        fastqc_zip_pre = FASTQC_RAW.out.zip
        fastqc_zip_post = FASTQC_CLEAN.out.zip
        trimmed = TRIMGALORE.out.reads
        clean = KNEADDATA.out.reads
        summary_plot = SUMMARY_KNEADDATA.out.qc_summary_plot
        versions = ch_versions_qc
}
