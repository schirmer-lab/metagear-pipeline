include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"
include { createExistingDirChannel; createExistingFileChannel } from "$projectDir/subworkflows/local/utils/existing_data"

include { ASSEMBLY } from "$projectDir/subworkflows/local/common/assembly"

include { GENE_CALL } from "$projectDir/subworkflows/local/common/gene_call"
include { PROTEIN_CALL } from "$projectDir/subworkflows/local/common/protein_call"
include { PROTEIN_ANNOTATION } from "$projectDir/subworkflows/local/common/protein_annotation"

include { BWA_INDEX } from "$projectDir/modules/nf-core/bwa/index"
include { ABUNDANCE } from "$projectDir/subworkflows/local/common/abundance"

include { CLUSTER_SEQUENCES as CLUSTER_GENES; CLUSTER_SEQUENCES as CLUSTER_PROTEINS } from "$projectDir/subworkflows/local/common/clustering"

workflow GENES_INIT {

    main:
        if ( params.input ) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

        amrfinder_db = Channel.fromPath("${params.amrfinder_db}", checkIfExists: true)

        INPUT_CHECK ( ch_input, "reads" )

    emit:
        validated_input = INPUT_CHECK.out.validated_input
        amrfinder_db
        versions = INPUT_CHECK.out.versions
}


workflow GENES {

    take:
        clean_reads // [meta, reads]
        amrfinder_db

    main:

        ch_versions = Channel.empty()

        ch_contigs = Channel.empty()

        if ( params.contigs_dir ) {
            ch_contigs = createExistingDirChannel ( params.contigs_dir, "*.contigs.fa.gz", ".contigs.fa", false )

        }else {
            // Assemble clean_reads into contigs
            ASSEMBLY ( clean_reads )
            ch_versions = ch_versions.mix( ASSEMBLY.out.versions.first() )

            ch_contigs = ASSEMBLY.out.contigs
        }

        ch_genes = Channel.empty()
        if ( params.genes_dir ) {
            ch_genes = createExistingDirChannel ( params.genes_dir, "*.all.genes.filtered.fasta", ".all.genes.filtered", { [ [id: it[0].id + '.all.genes', label: 'all.genes', src: it[0].id ], it[1] ] } )

        } else {
            // Call genes for all contigs
            ch_gene_call = ch_contigs.map { [ [id: it[0].id + '.all.genes', label: 'all.genes', src: it[0].id ], it[1] ] }

            GENE_CALL ( ch_gene_call )
            ch_versions = ch_versions.mix( GENE_CALL.out.versions.first() )

            ch_genes = GENE_CALL.out.genes
        }

        ch_representative_genes = Channel.empty()

        if ( params.representative_genes ) {
            ch_representative_genes = createExistingFileChannel ( params.representative_genes, { [ [id: "all.genes"], it ] } )

        } else {
            CLUSTER_GENES ( ch_genes, "mmseqs2", true )
            ch_versions =  ch_versions.mix(CLUSTER_GENES.out.versions)
            ch_representative_genes = CLUSTER_GENES.out.representative
        }

        ch_representative_proteins = Channel.empty()

        if ( params.representative_proteins ) {
            ch_representative_proteins = createExistingFileChannel ( params.representative_proteins, { [ [id: "all.genes"], it ] } )

        } else {
            PROTEIN_CALL ( ch_representative_genes )
            ch_versions =  ch_versions.mix(PROTEIN_CALL.out.versions)

            CLUSTER_PROTEINS ( PROTEIN_CALL.out.proteins, "mmseqs2", false )
            ch_representative_proteins = CLUSTER_PROTEINS.out.representative
        }

        if ( !params.representative_proteins_annotations ) {
            PROTEIN_ANNOTATION ( ch_representative_proteins, amrfinder_db )
            ch_versions =  ch_versions.mix(PROTEIN_ANNOTATION.out.versions)
        }

        ch_representative_genes_tpm = Channel.empty()
        ch_representative_genes_rpkm = Channel.empty()
        ch_representative_genes_count = Channel.empty()

        if ( params.representative_genes_tpm && params.representative_genes_rpkm && params.representative_genes_count ) {

            ch_representative_genes_tpm = createExistingFileChannel ( params.representative_genes_tpm, { [ [id: "all.genes_tpm"], it ] } )
            ch_representative_genes_rpkm = createExistingFileChannel ( params.representative_genes_rpkm, { [ [id: "all.genes_rpkm"], it ] } )
            ch_representative_genes_count = createExistingFileChannel ( params.representative_genes_count, { [ [id: "all.genes_count"], it ] } )

        } else {
            ch_abundance_input = clean_reads.combine( ch_representative_genes )
                                .map { meta_reads, reads, meta_genes, genes -> [ meta_reads + [label: meta_genes.id], reads, genes ] }

            ABUNDANCE ( ch_abundance_input, 'contig', file("$projectDir/assets/empty.txt") )
            ch_versions =  ch_versions.mix(ABUNDANCE.out.versions)

            ch_representative_genes_count = ABUNDANCE.out.count
            ch_representative_genes_rpkm = ABUNDANCE.out.rpkm
            ch_representative_genes_tpm = ABUNDANCE.out.tpm
        }

    emit:
        // TODO: implement emission of all other relevant channels
        versions = ch_versions
}
