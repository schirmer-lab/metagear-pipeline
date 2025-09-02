include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { ASSEMBLY } from "$projectDir/subworkflows/local/common/assembly"

include { VIRAL_DETECTION } from "$projectDir/subworkflows/local/virus/detection"
include { VIRAL_ANNOTATION } from "$projectDir/subworkflows/local/virus/annotation"

include {BWA_INDEX} from "$projectDir/modules/nf-core/bwa/index"

include { GENE_CALL } from "$projectDir/subworkflows/local/common/gene_call"
include { PROTEIN_CALL; PROTEIN_CALL as PLASMID_PROTEIN_CALL } from "$projectDir/subworkflows/local/common/protein_call"

include { PROTEIN_ANNOTATION } from "$projectDir/subworkflows/local/common/protein_annotation"

// include { ABUNDANCE as VOTU_ABUNDANCE; ABUNDANCE as GENE_ABUNDANCE; ABUNDANCE as PLASMID_ABUNDANCE; ABUNDANCE as PLASMID_GENE_ABUNDANCE } from "$projectDir/subworkflows/local/common/abundance"
include { ABUNDANCE } from "$projectDir/subworkflows/local/common/abundance"

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

        // Modify sequences to append "viral_genes|plasmid_genes"
        ch_gene_call = VIRAL_DETECTION.out.sequences.map { [ [id: it[0].id + '.' + it[0].label + '.genes', label: it[0].label, src: it[0].id ], it[1] ] }

        GENE_CALL ( ch_gene_call )

        // Build index for abundance estimation (only once)
        BWA_INDEX ( VIRAL_DETECTION.out.catalogs.concat(GENE_CALL.out.gene_catalog) )

        ch_abundance_input = reads.combine( BWA_INDEX.out.index )
                                .map { meta_reads, reads, meta_index, index -> [ [id: meta_reads.id, label: meta_index.id ], reads, index ] }

        ABUNDANCE ( ch_abundance_input )

        PROTEIN_CALL ( GENE_CALL.out.gene_catalog )

        PROTEIN_ANNOTATION ( PROTEIN_CALL.out.protein_catalog )

        ch_protein_catalog = PROTEIN_CALL.out.protein_catalog
                                    .filter { meta, _ -> meta.id == 'virus.proteins' }
                                    .map { meta, fa -> [ [id: meta.id, is_proteins: true], fa ] }


        VIRAL_ANNOTATION ( VIRAL_DETECTION.out.viral_catalog, ch_protein_catalog, virsorter2_db, dram_db, iphop_db, amrfinder_db )

        // summary channel version
        ch_versions = ASSEMBLY.out.versions
                        .mix(VIRAL_DETECTION.out.versions)
                        .mix(GENE_CALL.out.versions)
                        .mix(ABUNDANCE.out.versions)
                        .mix(PROTEIN_CALL.out.versions)
                        .mix(PROTEIN_ANNOTATION.out.versions)
                        .mix(VIRAL_ANNOTATION.out.versions)

    emit:

        versions = ch_versions

}
