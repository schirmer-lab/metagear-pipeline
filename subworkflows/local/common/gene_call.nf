/* --- Assembly and Gene Calling --- */

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { MEGAHIT } from "$projectDir/modules/local/megahit/main"
include { PRODIGAL } from "$projectDir/modules/nf-core/prodigal"
include { FILTER_PRODIGAL } from "$projectDir/modules/local/metagear/utils/filter_prodigal"

include { CDHIT_CDHITEST } from "$projectDir/modules/local/cdhit/cdhitest/main"

/* --- Initialization for standalone process --- */
workflow GENE_CALL_INIT {
    main:

        if ( params.input ) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

        INPUT_CHECK ( ch_input, "reads" )

    emit:
        validated_input = INPUT_CHECK.out.validated_input
        versions = INPUT_CHECK.out.versions

}

/* --- Main Workflow --- */
workflow GENE_CALL {

    take:
        ch_clean_reads // meta, reads

    main:

        MEGAHIT (
            ch_clean_reads.map { meta, fastq -> [ meta, fastq ] }
        )

        PRODIGAL ( MEGAHIT.out.contigs, "gff" )

        FILTER_PRODIGAL ( PRODIGAL.out.nucleotide_fasta )

        // FILTER_PRODIGAL.out.filtered_fasta

        // CDHIT_CDHITEST (  )

        // ch_derreplicated_genes = CDHIT_CDHITEST.out.fasta

        ch_versions = MEGAHIT.out.versions.first()
                        .mix(PRODIGAL.out.versions.first())
                        .mix(CDHIT_CDHITEST.out.versions)


    emit:
        contigs = MEGAHIT.out.contigs
        versions = ch_versions
}
