include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { ASSEMBLY } from "$projectDir/subworkflows/local/common/assembly"

include { VIRAL_DETECTION } from "$projectDir/subworkflows/local/virus/detection"
include { VIRAL_ANNOTATION } from "$projectDir/subworkflows/local/virus/annotation"

include { GENE_CALL } from "$projectDir/subworkflows/local/common/gene_call"
include { PROTEIN_CALL } from "$projectDir/subworkflows/local/common/protein_call"

include { PROTEIN_ANNOTATION } from "$projectDir/subworkflows/local/common/protein_annotation"

// include { VAMB_CONCATENATE_FASTA } from "$projectDir/modules/local/vamb/main"
// include { CLUSTER_SEQUENCES } from "$projectDir/subworkflows/local/common/clustering"
include { ABUNDANCE as VOTU_ABUNDANCE; ABUNDANCE as GENE_ABUNDANCE } from "$projectDir/subworkflows/local/common/abundance"


workflow VIRAL_ANALYSIS_INIT {

    main:
        if ( params.input ) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

        INPUT_CHECK ( ch_input, "reads" )

    emit:
        reads = INPUT_CHECK.out.validated_input
        versions = INPUT_CHECK.out.versions
}


workflow VIRAL_ANALYSIS {

    take:
        reads // [meta, reads(fast1, fas2)]

    main:

        ch_versions = Channel.empty()

        // ------ Assembly: Handle both missing contigs and full assembly ---- //
        ASSEMBLY ( reads )
        ch_versions = ch_versions.mix( ASSEMBLY.out.versions )

        genomad_db = Channel.fromPath("${params.genomad_db}", checkIfExists: true).first()
        checkv_db = Channel.fromPath("${params.checkv_db}", checkIfExists: true).first()

        virsorter2_db = Channel.fromPath("${params.virsorter2_db}", checkIfExists: true).first()
        dram_db = Channel.fromPath("${params.dram_db}", checkIfExists: true).first()
        iphop_db = Channel.fromPath("${params.iphop_db}", checkIfExists: true).first()
        amrfinder_db = Channel.fromPath("${params.amrfinder_db}", checkIfExists: true).first()

        VIRAL_DETECTION ( ASSEMBLY.out.contigs, genomad_db, checkv_db )

        VOTU_ABUNDANCE ( "votu_abundance", reads, VIRAL_DETECTION.out.viral_catalog )

        GENE_CALL ( "votu_gene_catalog", VIRAL_DETECTION.out.viral_catalog )

        GENE_ABUNDANCE ( "votu_gene_abundance", reads, GENE_CALL.out.gene_catalog )

        PROTEIN_CALL ( "votu_protein_catalog", GENE_CALL.out.gene_catalog )

        PROTEIN_ANNOTATION ( PROTEIN_CALL.out.protein_catalog )

        ch_protein_catalog = PROTEIN_CALL.out.protein_catalog.map { it -> [ [id: "votu_proteins", is_proteins: true], it[1] ] }

        VIRAL_ANNOTATION ( VIRAL_DETECTION.out.viral_catalog, ch_protein_catalog, virsorter2_db, dram_db, iphop_db, amrfinder_db )

        // summary channel version
        ch_versions = ASSEMBLY.out.versions
                        .mix(VIRAL_DETECTION.out.versions)
                        .mix(VOTU_ABUNDANCE.out.versions)
                        .mix(GENE_CALL.out.versions)
                        .mix(GENE_ABUNDANCE.out.versions)
                        .mix(PROTEIN_CALL.out.versions)
                        .mix(PROTEIN_ANNOTATION.out.versions)
                        .mix(VIRAL_ANNOTATION.out.versions)

    emit:

        versions = ch_versions

}
