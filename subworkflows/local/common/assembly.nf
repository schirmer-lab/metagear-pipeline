include { MEGAHIT } from "$projectDir/modules/local/megahit/main"

/* --- Main Workflow --- */
workflow ASSEMBLY {

    take:
        ch_clean_reads // [meta, reads] - reads that need assembly

    main:
        // Reuse of existing contigs is handled by the callers via params.contigs_dir;
        // this subworkflow is only reached when there are reads to assemble.
        MEGAHIT ( ch_clean_reads.map { meta, fastq -> [ meta, fastq ] } )

    emit:
        contigs = MEGAHIT.out.contigs
        versions = MEGAHIT.out.versions.first()
}
