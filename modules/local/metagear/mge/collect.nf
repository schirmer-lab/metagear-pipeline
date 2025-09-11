process COLLECT_TABLES {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/25.08.11/python_3.10.sif':
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    tuple val(meta), path(tables) // [ [id: virus|virus.filtered|plasmid|plasmid.filtered|dramv] path(files...) ]

    output:
    tuple val(meta), path("*.tsv"), emit: collected_tables
    tuple val(meta), path("*.txt"), emit: collected_ids, optional: true

    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [[ "${prefix}" == "virus" || "${prefix}" == "virus.filtered" ]]; then
        # Join virus summary tables (geNomad, CheckV)
        head -1 ${tables[0]} > ${prefix}_checkv_summary_taxa.tsv

        # Bulletproof approach for many files - use find with -exec or printf with xargs
        printf '%s\\n' ${tables} | xargs cat | grep -v virus_id >> ${prefix}_checkv_summary_taxa.tsv

        if [[ "${prefix}" == "virus.filtered" ]]; then
            grep -v virus_id ${prefix}_checkv_summary_taxa.tsv | cut -f2 > ${prefix}_ids.tmp
            grep '|' ${prefix}_ids.tmp > ${prefix}_provirus_ids.txt
            grep -v '|' ${prefix}_ids.tmp > ${prefix}_virus_ids.txt
        fi
    fi

    if [[ "${prefix}" == "plasmid" ||  "${prefix}" == "plasmid.filtered" ]]; then
        # Join plasmid tables (geNomad)
        head -1 ${tables[0]} > ${prefix}_summary.tsv
        printf '%s\\n' ${tables} | xargs cat | grep -v seq_name >> ${prefix}_summary.tsv

        if [[ "${prefix}" == "plasmid.filtered" ]]; then
            grep -v seq_name ${prefix}_summary.tsv | cut -f1 > ${prefix}_ids.txt
        fi
    fi

    if [[ "${prefix}" == "amg" ]]; then
        # Join AMG summary tables (DRAM-V)
        head -1 ${tables[0]} > ${prefix}_summary.tsv
        printf '%s\\n' ${tables} | xargs cat | grep -v gene >> ${prefix}_summary.tsv
    fi

    if [[ "${prefix}" == "host.genus" || "${prefix}" == "host.genome" ]]; then
        # Join hosts summary tables (iPhOP)
        head -1 ${tables[0]} > ${prefix}_summary.tsv
        printf '%s\\n' ${tables} | xargs cat | grep -v Virus >> ${prefix}_summary.tsv
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version | head -n1 | cut -d' ' -f4))
    END_VERSIONS
    """
}
