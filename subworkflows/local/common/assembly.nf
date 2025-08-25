include { MEGAHIT } from "$projectDir/modules/local/megahit/main"

/* --- Main Workflow --- */
workflow ASSEMBLY {

    take:
        ch_clean_reads // [meta, reads] - reads that need assembly

    main:
        ch_versions = Channel.empty()

        // Check if existing contigs are provided via params
        ch_existing_contigs = false

        if ( params.contigs ) {
            def contigs_path = file(params.contigs)

            if ( contigs_path.isDirectory() ) {
                ch_existing_contigs = Channel.fromPath("${params.contigs}/*.contigs.fa.gz")
                    .map { file ->
                        def sample_id = file.baseName.replaceAll(/\.contigs\.fa$/, '')
                        [ [id: sample_id], file ]
                    }
            } else {
                exit 1, "The provided contigs path is not a directory: ${contigs_path}"
            }
        }

        // If no existing contigs provided, assemble all reads
        if ( !ch_existing_contigs ) {
            MEGAHIT ( ch_clean_reads.map { meta, fastq -> [ meta, fastq ] } )
            ch_all_contigs = MEGAHIT.out.contigs
            ch_versions = MEGAHIT.out.versions.first()
        } else {
            // Determine which reads need assembly (missing contigs)
            ch_joined = ch_clean_reads.join(ch_existing_contigs, remainder: true)

            ch_missing_contigs = ch_joined.filter { meta, reads, contig -> contig == null }
                                         .map { meta, reads, contig -> [meta, reads] }

            ch_with_contigs = ch_joined.filter { meta, reads, contig -> contig != null }
                                      .map { meta, reads, contig -> [meta, contig] }

            // Assemble only missing contigs (if any)
            MEGAHIT ( ch_missing_contigs.map { meta, fastq -> [ meta, fastq ] } )
            ch_assembled_contigs = MEGAHIT.out.contigs
            ch_versions = MEGAHIT.out.versions.first()

            // Combine existing contigs with newly assembled ones
            ch_all_contigs = ch_with_contigs.mix(ch_assembled_contigs)
        }

    emit:
        contigs = ch_all_contigs
        versions = ch_versions
}
