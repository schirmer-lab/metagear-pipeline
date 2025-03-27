/*--- Check input samplesheet and get read channels ---*/

include { SAMPLESHEET_CHECK; RENAME_FILES; } from "$projectDir/modules/local/metagear/samplesheet_check"

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
                    // reorder the `reads` array based on the `indexes` array, where each element in `indexes` corresponds to the desired position of the element in `reads`
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
        validated_input  // channel: [ val(meta), [ etc ] ]
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

    def fastq_meta = []
    if (input_type == "contig") {
        fastq_meta = [ meta, file(row.contig) ]
    }
    else if ( input_type == "blast_seqs" ) {
        fastq_meta = [ meta, file(row.query_sequence), file(row.search_database) ]
    }else{
        if (row.fastq_2?.trim()) {
            fastq_meta = [ meta, file(row.fastq_1), file(row.fastq_2) ]
        }else{
            fastq_meta = [ meta, file(row.fastq_1) ]
        }
    }

    return fastq_meta
}
