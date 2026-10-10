/* --- Perform quality control on reads --- */

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { FASTQC as FASTQC_RAW; FASTQC as FASTQC_CLEAN} from "$projectDir/modules/nf-core/fastqc/main"
include { TRIMGALORE } from "$projectDir/modules/nf-core/trimgalore/main"
include { KNEADDATA; PARSE_KNEADDATA; SUMMARY_KNEADDATA } from "$projectDir/modules/local/kneaddata/main"

// Per-assay name first, so qc_dna and qc_rna can both stay configured; kneaddata_refdb
// remains the fallback for installations that already set it.
def qc_reference_databases(String assay) {
    def configured = params["qc_${assay}_refdb"] ?: params.kneaddata_refdb
    def paths = (configured instanceof List ? configured : [configured]).findAll { it }

    if ( !paths ) {
        error "qc_${assay}: no reference database configured. Set qc_${assay}_refdb (or kneaddata_refdb) to one or more bowtie2 index prefixes."
    }

    return paths.collect { file(it, checkIfExists: true) }
}

/* --- Main Workflow --- */
workflow QUALITY_CONTROL_INIT {
    take:
        assay // str: dna | rna

    main:

        if ( params.input ) { ch_input = file(params.input) } else { error 'Input samplesheet not specified!' }

        INPUT_CHECK ( ch_input, "reads" )

        kneaddata_refdb = Channel.value( qc_reference_databases(assay) )

    emit:
        validated_input = INPUT_CHECK.out.validated_input
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

        KNEADDATA ( TRIMGALORE.out.reads, kneaddata_refdb )

        PARSE_KNEADDATA ( KNEADDATA.out.kneaddata_log )

        PARSE_KNEADDATA.out.read_counts
            .map { meta, read_counts -> read_counts }
            .collect()
            .set { all_read_counts }

        SUMMARY_KNEADDATA ( all_read_counts )

        // Append "_clean" to IDs for post cleaning fastqc
        ch_fastqc_clean = KNEADDATA.out.reads
                            .map { meta, reads -> [[id: meta.id + "_clean"] , reads] }

        FASTQC_CLEAN ( ch_fastqc_clean )

        ch_versions_qc = FASTQC_RAW.out.versions
                            .mix(TRIMGALORE.out.versions)
                            .mix(KNEADDATA.out.versions)
                            .mix(FASTQC_CLEAN.out.versions)
                            .mix(PARSE_KNEADDATA.out.versions)
                            .mix(SUMMARY_KNEADDATA.out.versions)

    emit:
        fastqc_zip_pre = FASTQC_RAW.out.zip
        fastqc_zip_post = FASTQC_CLEAN.out.zip
        trimmed = TRIMGALORE.out.reads
        clean = KNEADDATA.out.reads
        read_counts = SUMMARY_KNEADDATA.out.read_counts
        summary_sections = SUMMARY_KNEADDATA.out.mqc
        versions = ch_versions_qc
}
