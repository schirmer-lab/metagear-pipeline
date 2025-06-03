import java.nio.file.Files
import java.nio.file.StandardCopyOption

include { GENOMAD_ENDTOEND as GENOMAD_PASS1; GENOMAD_ENDTOEND as GENOMAD_PASS2 } from "$projectDir/modules/nf-core/genomad/endtoend"
include { CHECKV_ENDTOEND as CHECKV_PASS1; CHECKV_ENDTOEND as CHECKV_PASS2} from "$projectDir/modules/nf-core/checkv/endtoend"

include { CHECKV_ADAPT_OUTPUT } from "$projectDir/modules/local/checkv/adapt"
include { CREATE_TABLES } from "$projectDir/modules/local/mvip/create_tables"

workflow VIRAL_DETECTION {

    take:
        contigs // [ id, fasta] -> [ [id: vmx.bins], [fasta.gz] ]
        genomad_db
        checkv_db

    main:

        // 1. Initial geNomad pass
        GENOMAD_PASS1 ( contigs, genomad_db )
        ch_pass1_viruses = GENOMAD_PASS1.out.virus_fasta.filter { file -> file.size() > 0 }

        // 2. CheckV on viruses from geNomad
        CHECKV_PASS1 (ch_pass1_viruses, checkv_db)
        ch_pass1_proviruses = CHECKV_PASS1.out.proviruses.filter { file -> file.size() > 0 }

        // 3. Second geNomad pass on trimmed provirus
        GENOMAD_PASS2( ch_pass1_proviruses, genomad_db )
        ch_pass2_viruses = GENOMAD_PASS2.out.virus_fasta.filter { file -> file.size() > 0 }

        // 4. Final CheckV on non-empty geNomad2-detected viral contigs
        checkv2 = CHECKV_PASS2 ( ch_pass2_viruses, checkv_db )

        // 5. Collect all summary files (2 or 4) for each sample
        ch_pass1_virus_summary = GENOMAD_PASS1.out.virus_summary.map { [ it[0], [ it[1], "virus" ] ] }
        ch_checkv1_summary = CHECKV_PASS1.out.quality_summary.map { [ it[0], [ it[1], "virus" ] ] }
        ch_pass2_virus_summary = GENOMAD_PASS2.out.virus_summary.map { [ it[0], [ it[1], "provirus" ] ] }
        ch_checkv2_summary = CHECKV_PASS2.out.quality_summary.map { [ it[0], [ it[1], "provirus" ] ] }


        ch_all_summaries = ch_pass1_virus_summary
                            .mix(ch_checkv1_summary)
                            .mix(ch_pass2_virus_summary)
                            .mix(ch_checkv2_summary)
                            .groupTuple( by: 0 )

        ch_all_summaries.view()

        def ictv_taxonomy = Channel.fromPath("$projectDir/assets/metagear/ICTV_Taxonomy_List.tsv", checkIfExists: true).first()
        CREATE_TABLES ( ch_all_summaries, ictv_taxonomy )
        // ch_grouped = ch_all_summaries
        //     .groupTuple {  }   // group by sample_id
        //     .map { sample, items -> tuple(sample, items.collect { it[1] }) }

        // ch_grouped.view()
        // Feed grouped files into CREATE_TABLES
        // CREATE_TABLES(ch_grouped)


    emit:
        //TODO:...
        versions = []
}
