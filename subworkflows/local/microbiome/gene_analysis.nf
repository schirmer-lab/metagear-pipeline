include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { ASSEMBLY } from "$projectDir/subworkflows/local/common/assembly"

include { GENE_CALL } from "$projectDir/subworkflows/local/common/gene_call"

include {BWA_INDEX} from "$projectDir/modules/nf-core/bwa/index"

include { PROTEIN_CALL } from "$projectDir/subworkflows/local/common/protein_call"

include { ABUNDANCE } from "$projectDir/subworkflows/local/common/abundance"

include { PROTEIN_ANNOTATION } from "$projectDir/subworkflows/local/common/protein_annotation"

include { MSP } from "$projectDir/subworkflows/local/pangenome/msp"

workflow GENE_ANALYSIS_INIT {

    main:
        if ( params.input ) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

        gtdb_tk_db = Channel.fromPath("${params.gtdb_tk_db}", checkIfExists: true)
        amrfinder_db = Channel.fromPath("${params.amrfinder_db}", checkIfExists: true)

        metaphlan_db = Channel.empty()

        if ( params.metaphlan_profiles ) {
            metaphlan_profiles = Channel.fromPath("${params.metaphlan_profiles}", checkIfExists: true)
        } else {
            metaphlan_db = Channel.fromPath("${params.metaphlan_db}", checkIfExists: true).first()
            metaphlan_profiles = false
        }

        INPUT_CHECK ( ch_input, "reads" )

    emit:
        validated_input = INPUT_CHECK.out.validated_input
        metaphlan_profiles
        gtdb_tk_db
        metaphlan_db
        amrfinder_db
        versions = INPUT_CHECK.out.versions
}


workflow GENE_ANALYSIS {

    take:
        clean_reads // [meta, reads]
        metaphlan_profiles
        gtdb_tk_db
        amrfinder_db

    main:

        ASSEMBLY ( clean_reads )

        ch_gene_call = ASSEMBLY.out.contigs.map { [ [id: it[0].id + '.all.genes', label: 'all', src: it[0].id ], it[1] ] }

        GENE_CALL ( ch_gene_call )

        BWA_INDEX ( GENE_CALL.out.gene_catalog )

        ch_abundance_input = clean_reads.combine( BWA_INDEX.out.index )
                                .map { meta_reads, reads, meta_index, index -> [ [id: meta_reads.id, label: meta_index.id ], reads, index ] }

        ABUNDANCE ( ch_abundance_input )

        PROTEIN_CALL ( GENE_CALL.out.gene_catalog )

        PROTEIN_ANNOTATION ( PROTEIN_CALL.out.protein_catalog, amrfinder_db )

        MSP ( GENE_CALL.out.gene_catalog, ABUNDANCE.out.count, ABUNDANCE.out.rpkm, gtdb_tk_db, metaphlan_profiles )

        // summary channel version
        ch_versions = ASSEMBLY.out.versions
                        .mix(GENE_CALL.out.versions)
                        .mix(ABUNDANCE.out.versions)
                        .mix(PROTEIN_CALL.out.versions)
                        .mix(PROTEIN_ANNOTATION.out.versions)
                        // .mix(MSP.out.versions)


    emit:
        // TODO: implement emission of all relevant channels
        versions = ch_versions
}
