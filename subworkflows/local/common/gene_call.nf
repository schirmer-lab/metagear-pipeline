/* --- Assembly and Gene Calling --- */

include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { PRODIGAL } from "$projectDir/modules/nf-core/prodigal"
include { FILTER_PRODIGAL } from "$projectDir/modules/local/metagear/utils/filter_prodigal"

include { VAMB_CONCATENATE_FASTA } from "$projectDir/modules/local/vamb/main"

// include { CDHIT_CDHITEST } from "$projectDir/modules/local/cdhit/cdhitest/main"
include { MMSEQS_EASY_CLUSTER } from "$projectDir/modules/local/mmseqs/easy_cluster/main"

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
        
        // Build gene catalog from scratch
        PRODIGAL ( sequences, "gff" )

        FILTER_PRODIGAL ( PRODIGAL.out.nucleotide_fasta )

        ch_merged_genes = FILTER_PRODIGAL.out.filtered_fasta.map{ [[ id: it[0].label + '.genes' ], it[1] ] }
                .groupTuple(by: 0)

        VAMB_CONCATENATE_FASTA ( ch_merged_genes )

        // CDHIT_CDHITEST ( VAMB_CONCATENATE_FASTA.out.catalog )
        MMSEQS_EASY_CLUSTER ( VAMB_CONCATENATE_FASTA.out.catalog )

        ch_gene_catalog = MMSEQS_EASY_CLUSTER.out.representatives

        ch_versions = PRODIGAL.out.versions.first()
                        .mix(FILTER_PRODIGAL.out.versions.first())
                        .mix(VAMB_CONCATENATE_FASTA.out.versions)
                        .mix(MMSEQS_EASY_CLUSTER.out.versions)
    

    emit:
        gene_catalog = ch_gene_catalog
        versions = ch_versions
}
