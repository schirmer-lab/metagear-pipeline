include { MEGAHIT } from "$projectDir/modules/local/megahit/main"

/* --- Main Workflow --- */
workflow ASSEMBLY {

    take:
        ch_clean_reads // meta, reads

    main:

        MEGAHIT (
            ch_clean_reads.map { meta, fastq -> [ meta, fastq ] }
        )

        ch_versions = MEGAHIT.out.versions.first()

    emit:
        contigs = MEGAHIT.out.contigs
        versions = ch_versions
}
