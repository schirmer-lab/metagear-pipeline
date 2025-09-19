/* --- Assembly and Gene Calling --- */

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { PRODIGAL } from "$projectDir/modules/nf-core/prodigal"
include { FILTER_PRODIGAL } from "$projectDir/modules/local/metagear/utils/filter_prodigal"

include { EXTRACT_GENES } from "$projectDir/modules/local/metagear/mge/extract_genes"

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
        sequences // tuple (meta, reads) -> [ [id: sample1_*plasmid|votu*, label: plasmid_genes|votu_genes|genes ], fasta ]

    main:

        // sequences.map { meta, _ -> meta.id = meta.label ? meta.id + '.' + meta.label : meta.id + '.genes' }

        // Build gene catalog from scratch
        PRODIGAL ( sequences, "gff" )

        FILTER_PRODIGAL ( PRODIGAL.out.nucleotide_fasta )

        ch_versions = PRODIGAL.out.versions.first()
                        .mix(FILTER_PRODIGAL.out.versions.first())


    emit:
        genes = FILTER_PRODIGAL.out.filtered_fasta
        versions = ch_versions
}

workflow VIRAL_GENE_CALL {

    take:
        viral_sequences // tuple (meta, reads) -> [ [id: sample1_*plasmid|votu*, label: plasmid_genes|votu_genes|genes ], fasta ]

    main:
        ch_versions = Channel.empty()

        EXTRACT_GENES ( viral_sequences )
        ch_versions =  ch_versions.mix( EXTRACT_GENES.out.versions)

    emit:
        versions = ch_versions
        genes = EXTRACT_GENES.out.extracted_genes
}