/*--- Check input samplesheet and get read channels ---*/

import groovy.transform.Field

include { SAMPLESHEET_CHECK; RENAME_FILES; } from "$projectDir/modules/local/metagear/samplesheet_check"

// Allowed biome values (SemiBin2 pretrained environments + 'global' fallback).
// Keep in sync with assets/schema_input.json. Validated in Groovy (not in
// bin/input_validator.py) so adding/removing biomes doesn't invalidate the
// per-task bin/ hash that Nextflow uses for caching.
//
// @Field is required: without it, `def BIOMES` only binds in the script's
// main flow, and the script-scoped `create_input_channel` function below
// cannot see it (Nextflow + Groovy script-scope quirk — surfaces at runtime
// as "No such property: BIOMES" when the .map closure invokes the function).
@Field
final List BIOMES = [
    'human_gut', 'human_oral', 'mouse_gut', 'dog_gut', 'cat_gut',
    'ocean', 'soil', 'built_environment', 'wastewater', 'chicken_caecum',
    'global'
]

workflow INPUT_CHECK {
    take:
        samplesheet // file: /path/to/samplesheet.csv
        input_type // str: [reads, contig, contig_reads, blast_seqs]

    main:

        csv_channel = SAMPLESHEET_CHECK ( samplesheet, input_type ).csv.splitCsv ( header:true, sep:',' )

        csv_channel.map { create_input_channel(it, input_type) }

                    .filter{ !it[0].id.startsWith("#") } // Filter out lines starting with '#'
                    .flatMap { it ->
                            def tuples = [ [it[0], 0, it[1]] ]
                            if (it[2]){
                                tuples.add( [it[0], 1, it[2]])
                            }
                            if (it[3]){
                                tuples.add( [it[0], 2, it[3]])
                            }
                            return tuples
                    }
                    .set { input_tuples }

        RENAME_FILES ( input_tuples )

        RENAME_FILES.out.renamed_files
            .groupTuple(by: 0)
            .map{ it ->
                    def indexes = it[1]
                    def reads = it[2]
                    ordered_reads = []

                    for (int i = 0; i < indexes.size(); i++) {
                        for (int j = 0; j < reads.size(); j++) {
                            if (indexes[j] == i) {
                                ordered_reads.add(reads[j])
                            }
                        }
                    }

                    return [ it[0], ordered_reads ]
            }
            .set { validated_input }

    emit:
        validated_input = validated_input  // channel: [ val(meta), [ etc ] ]
        versions = SAMPLESHEET_CHECK.out.versions.mix(RENAME_FILES.out.versions) // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_input_channel(LinkedHashMap row, String input_type) {
    // create meta map
    def meta = [:]

    if (input_type == "blast_seqs") {
        meta.id = row.analysis
    }else{
        meta.id = row.sample
    }

    if (input_type == "grouped_reads") {
        meta.group = row.group
        meta.tag = row.tag
    }

    // The samplesheet's optional `biome` column is validated here (against the
    // SemiBin2 environment list) but intentionally NOT copied into meta.
    // Adding a key to meta would change every downstream task's cache hash for
    // existing workflows that don't need biome (genes, virus,
    // microbial_profiles). The future bacterial_binning subworkflow will read
    // biome from the CSV directly in its _INIT and join it onto its reads
    // channel locally, scoping the meta change to just its own processes.
    def biome_value = row.biome?.trim()
    if (biome_value && !BIOMES.contains(biome_value)) {
        exit 1, "Invalid biome '${biome_value}' for sample '${meta.id}'. Allowed values: ${BIOMES.join(', ')} (or leave blank to default to 'global')."
    }

    def fastq_meta = []

    if (input_type == "contig") {
        fastq_meta = [ meta, file(row.contig) ]
    }
    else if ( input_type == "blast_seqs" ) {
        fastq_meta = [ meta, file(row.query_sequence), file(row.search_database) ]
    }else if ( input_type == "contig_reads" ) {
        fastq_meta = [ meta, file(row.contig), file(row.fastq_1), file(row.fastq_2) ]
    }else{
        if (row.fastq_2?.trim()) {
            fastq_meta = [ meta, file(row.fastq_1), file(row.fastq_2) ]
        }else{
            meta.single_end = true
            fastq_meta = [ meta, file(row.fastq_1) ]
        }
    }

    return fastq_meta
}
