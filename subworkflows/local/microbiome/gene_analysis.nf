include { INPUT_CHECK } from "$projectDir/subworkflows/local/common/input_check"

include { ASSEMBLY } from "$projectDir/subworkflows/local/common/assembly"

include { GENE_CALL } from "$projectDir/subworkflows/local/common/gene_call"
include { PROTEIN_CALL } from "$projectDir/subworkflows/local/common/protein_call"
include { PROTEIN_ANNOTATION } from "$projectDir/subworkflows/local/common/protein_annotation"

include { BWA_INDEX } from "$projectDir/modules/nf-core/bwa/index"
include { ABUNDANCE } from "$projectDir/subworkflows/local/common/abundance"

include { CLUSTER_SEQUENCES as CLUSTER_GENES; CLUSTER_SEQUENCES as CLUSTER_PROTEINS } from "$projectDir/subworkflows/local/common/clustering"

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

        ch_versions = Channel.empty()

        ASSEMBLY ( clean_reads )
        ch_versions =  ch_versions.mix(ASSEMBLY.out.versions)

        ch_gene_call = ASSEMBLY.out.contigs.map { [ [id: it[0].id, label: 'all.genes', src: it[0].id ], it[1] ] }

        GENE_CALL ( ch_gene_call )
        ch_versions =  ch_versions.mix(GENE_CALL.out.versions)

        CLUSTER_GENES ( GENE_CALL.out.genes, "mmseqs2", true )
        ch_versions =  ch_versions.mix(CLUSTER_GENES.out.versions)

        PROTEIN_CALL ( CLUSTER_GENES.out.representative  )
        ch_versions =  ch_versions.mix(PROTEIN_CALL.out.versions)

        CLUSTER_PROTEINS ( PROTEIN_CALL.out.proteins, "mmseqs2", false )

        PROTEIN_ANNOTATION ( CLUSTER_PROTEINS.out.representative, amrfinder_db )
        ch_versions =  ch_versions.mix(PROTEIN_ANNOTATION.out.versions)

        BWA_INDEX ( CLUSTER_GENES.out.representative )
        ch_versions =  ch_versions.mix(BWA_INDEX.out.versions)

        ch_abundance_input = clean_reads.combine( BWA_INDEX.out.index )
                                .map { meta_reads, reads, meta_index, index -> [ [id: meta_reads.id, label: meta_index.id ], reads, index ] }

        ABUNDANCE ( ch_abundance_input )
        ch_versions =  ch_versions.mix(ABUNDANCE.out.versions)

        MSP ( CLUSTER_GENES.out.representative, ABUNDANCE.out.count, ABUNDANCE.out.rpkm, gtdb_tk_db, metaphlan_profiles )
        // ch_versions =  ch_versions.mix(MSP.out.versions) //TODO: Needs fixing, not working properly

    emit:
        // TODO: implement emission of all other relevant channels
        versions = ch_versions
}
