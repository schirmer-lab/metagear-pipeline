import java.nio.file.Files
import java.nio.file.StandardCopyOption

include { GENOMAD_ENDTOEND as GENOMAD_PASS1; GENOMAD_ENDTOEND as GENOMAD_PASS2 } from "$projectDir/modules/nf-core/genomad/endtoend"
include { CHECKV_ENDTOEND as CHECKV_PASS1; CHECKV_ENDTOEND as CHECKV_PASS2} from "$projectDir/modules/nf-core/checkv/endtoend"

include { CHECKV_ADAPT_OUTPUT } from "$projectDir/modules/local/checkv/adapt"
include { MERGE_VIRUS_TABLES } from "$projectDir/modules/local/mvip/create_tables"

include { SEQTK_SUBSEQ } from "$projectDir/modules/local/seqtk/subseq/main"

workflow VIRAL_DETECTION {

    take:
        contigs // [ id, fasta] -> [ [id: vmx.bins], [fasta.gz] ]
        genomad_db
        checkv_db

    main:

        // 1. Initial geNomad pass
        GENOMAD_PASS1 ( contigs, genomad_db )
        ch_pass1_viruses = GENOMAD_PASS1.out.virus_fasta.filter { meta, fna -> fna.toFile().length() > 0 }

        // 2. CheckV on viruses from geNomad
        CHECKV_PASS1 (ch_pass1_viruses, checkv_db)
        ch_pass1_proviruses = CHECKV_PASS1.out.proviruses.filter { meta, fna -> fna.toFile().length() > 0 }
    
        // 3. Second geNomad pass on trimmed provirus
        GENOMAD_PASS2( ch_pass1_proviruses, genomad_db )
        // ch_pass2_viruses = GENOMAD_PASS2.out.virus_fasta.filter { meta, fna -> fna.toFile().length() > 0 }
        ch_pass2_viruses = GENOMAD_PASS2.out.virus_fasta.filter { meta, gz ->
                                new java.util.zip.GZIPInputStream(gz.toFile().newInputStream())
                                    .withCloseable { it.read() != -1 }
                            }

        // 4. Final CheckV on non-empty geNomad2-detected viral contigs
        checkv2 = CHECKV_PASS2 ( ch_pass2_viruses, checkv_db )

        // 5. Adapt CheckV and geNomad channel outputs to allow merging
        ch_all_summaries = GENOMAD_PASS1.out.virus_summary.map { [ it[0], [ it[1], "virus" ] ] }
                            .mix(CHECKV_PASS1.out.quality_summary.map { [ it[0], [ it[1], "virus" ] ] })
                            .mix(GENOMAD_PASS2.out.virus_summary.map { [ it[0], [ it[1], "provirus" ] ] })
                            .mix(CHECKV_PASS2.out.quality_summary.map { [ it[0], [ it[1], "provirus" ] ] })
                            .groupTuple( by: 0 )

        def ictv_taxonomy = Channel.fromPath("$projectDir/assets/metagear/ICTV_Taxonomy_List.tsv", checkIfExists: true).first()
        MERGE_VIRUS_TABLES ( ch_all_summaries, ictv_taxonomy )

        joint = GENOMAD_PASS1.out.virus_fasta
            .mix ( GENOMAD_PASS2.out.virus_fasta )
            .groupTuple( by: 0 )
            .join ( MERGE_VIRUS_TABLES.out.sequence_ids, by: 0 )

        SEQTK_SUBSEQ ( joint )

    ch_versions = GENOMAD_PASS1.out.versions.first()
                    .mix(CHECKV_PASS1.out.versions.first())
                    .mix(SEQTK_SUBSEQ.out.versions.first())
                    // .mix(MERGE_VIRUS_TABLES.out.versions.first()) //TODO: versions is not being generated correctly, skipping for now


    emit:
        merged_tables   = MERGE_VIRUS_TABLES.out.merged_tables
        filtered_tables = MERGE_VIRUS_TABLES.out.filtered_tables
        sequences       = SEQTK_SUBSEQ.out.sequences
        versions        = ch_versions
}
