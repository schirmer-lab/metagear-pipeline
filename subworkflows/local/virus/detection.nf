import java.nio.file.Files
import java.nio.file.StandardCopyOption

include { MMSEQS_EASY_CLUSTER } from "$projectDir/modules/local/mmseqs/easy_cluster/main"

include { GENOMAD_ENDTOEND as GENOMAD_PASS1; GENOMAD_ENDTOEND as GENOMAD_PASS2 } from "$projectDir/modules/nf-core/genomad/endtoend"
include { CHECKV_ENDTOEND as CHECKV_PASS1; CHECKV_ENDTOEND as CHECKV_PASS2} from "$projectDir/modules/nf-core/checkv/endtoend"

include { CHECKV_ADAPT_OUTPUT } from "$projectDir/modules/local/checkv/adapt"
include { MERGE_TABLES } from "$projectDir/modules/local/mvip/create_tables"


// include { SEQTK_SUBSEQ as CONCAT_VIRUS; SEQTK_SUBSEQ as CONCAT_PLASMIDS} from "$projectDir/modules/local/seqtk/subseq/main"
include { SEQTK_SUBSEQ as EXTRACT_SEQUENCES } from "$projectDir/modules/local/seqtk/subseq/main"
include { VAMB_CONCATENATE_FASTA } from "$projectDir/modules/local/vamb/main"

include { COLLECT_TABLES } from "$projectDir/modules/local/metagear/mge/summarize"

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
                            .groupTuple( by: 0, size: 4, remainder: true )
                            .join(GENOMAD_PASS1.out.plasmid_summary)

        def ictv_taxonomy = Channel.fromPath("$projectDir/assets/metagear/ICTV_Taxonomy_List.tsv", checkIfExists: true).first()

        // Merge summary tables from geNomad and CheckV, apply taxonomy. Filter plasmids using FDR
        MERGE_TABLES ( ch_all_summaries, ictv_taxonomy )

        virus_to_keep = GENOMAD_PASS1.out.virus_fasta
                        .mix ( GENOMAD_PASS2.out.virus_fasta )
                        .groupTuple( by: 0, size: 2, remainder: true )
                        .join ( MERGE_TABLES.out.sequence_ids, by: 0 )
                        .map { [ [id: it[0].id, label: 'virus'], it[1], it[2] ] }

        plasmids_to_keep = GENOMAD_PASS1.out.plasmid_fasta
                            .join ( MERGE_TABLES.out.plasmid_sequence_ids, by: 0 )
                            .map { [ [id: it[0].id, label: 'plasmid'], it[1], it[2] ] }

        // Extract sequences to keep after filtering
        EXTRACT_SEQUENCES ( virus_to_keep.concat ( plasmids_to_keep ) )

        // Concatenate all viral sequences for clustering (derreplicated catalogs)
        ch_catalog_input = EXTRACT_SEQUENCES.out.sequences
                                .map { meta, fasta -> [ [id: meta.label], fasta ] }   // normalize meta
                                .groupTuple(by: 0) // collect all for the same id
                                .map { meta, paths -> [ meta, paths.sort { it.toString() } ] }

        VAMB_CONCATENATE_FASTA ( ch_catalog_input )

        MMSEQS_EASY_CLUSTER ( VAMB_CONCATENATE_FASTA.out.catalog )

        ch_versions = GENOMAD_PASS1.out.versions.first()
                        .mix(CHECKV_PASS1.out.versions.first())
                        .mix(EXTRACT_SEQUENCES.out.versions)
                    // .mix(MERGE_TABLES.out.versions.first()) //TODO: versions is not being generated correctly, skipping for now


    emit:
        viral_sequences =  EXTRACT_SEQUENCES.out.sequences.filter { meta, _ -> meta.label == 'virus' }
        viral_catalog = MMSEQS_EASY_CLUSTER.out.representatives.filter { meta, _ -> meta.id == 'virus' }
        plasmid_sequences = EXTRACT_SEQUENCES.out.sequences.filter { meta, _ -> meta.label == 'plasmid' }
        plasmid_catalog = MMSEQS_EASY_CLUSTER.out.representatives.filter { meta, _ -> meta.id == 'plasmid' }
        sequences = EXTRACT_SEQUENCES.out.sequences
        catalogs = MMSEQS_EASY_CLUSTER.out.representatives
        versions = ch_versions
}
