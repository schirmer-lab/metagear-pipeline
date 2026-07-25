process BUILD_MAG_CATALOG {
    tag "cohort"
    label 'process_low'

    conda "bioconda::seqkit=2.8.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.0--h9ee0642_0':
        'biocontainers/seqkit:2.8.0--h9ee0642_0' }"

    // Concatenates the cluster representatives from DREP into a single cohort
    // catalog FASTA and emits a contig→genome lookup TSV that
    // `coverm genome --genome-definition` consumes downstream.
    //
    // Across samples, raw assembly contig names collide (every MEGAHIT run
    // emits k141_0, k141_1, …), so each header is prefixed with the MAG name:
    //
    //     >k141_42   (in SAMPLE-0.binette_bin1.fa)
    //   →
    //     >SAMPLE-0.binette_bin1__k141_42
    //
    // Output TSVs:
    //   contig_to_mag.tsv     <mag>\t<renamed_contig>     (coverm input)
    //   mag_contig_lengths.tsv <renamed_contig>\t<length>  (sanity / dl-stats)

    input:
    path(reps, stageAs: 'reps/*')

    output:
    // Plain paths — Nextflow's parser refuses map literals inside val() at
    // module level ("No such variable: id"). Consumers wrap with meta as
    // needed in the caller subworkflow.
    path('mag_catalog.fa.gz'),      emit: catalog
    path('contig_to_mag.tsv'),      emit: contig_to_mag
    path('mag_contig_lengths.tsv'), emit: contig_lengths
    path 'versions.yml',            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    : > mag_catalog.fa
    : > metadata.tsv

    shopt -s nullglob
    for f in reps/*.fa reps/*.fasta reps/*.fna; do
        [ -f "\$f" ] || continue
        bn=\$(basename "\$f")
        # Strip the FASTA suffix to derive the MAG name (matches dRep/GTDB-Tk
        # cluster IDs which also drop the suffix).
        mag="\${bn%.fa}"
        mag="\${mag%.fasta}"
        mag="\${mag%.fna}"

        # 1. Rename headers and append to the cohort catalog.
        seqkit replace -p '^(\\S+)' -r "\${mag}__\\\$1" "\$f" >> mag_catalog.fa

        # 2. Pull contig name + length from the ORIGINAL file (fx2tab is
        #    cheaper than re-reading the catalog at the end). Prefix the
        #    contig name with mag__ here to match the renamed catalog.
        seqkit fx2tab -n -l "\$f" | awk -v OFS='\\t' -v mag="\$mag" '
            { print mag, mag "__" \$1, \$2 }
        ' >> metadata.tsv
    done

    # Split the metadata into the two consumer-ready TSVs.
    awk -F'\\t' 'BEGIN{OFS="\\t"} { print \$1, \$2 }' metadata.tsv > contig_to_mag.tsv
    awk -F'\\t' 'BEGIN{OFS="\\t"} { print \$2, \$3 }' metadata.tsv > mag_contig_lengths.tsv

    gzip mag_catalog.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """
}
