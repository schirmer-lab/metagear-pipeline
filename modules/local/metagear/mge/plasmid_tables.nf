process MERGE_PLASMID_TABLES {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/25.08.11/python_3.10.sif':
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    tuple val(meta), path(plasmid_tables)

    output:
    tuple val(meta), path("*_plasmid_summary.tsv"), emit: summary_plasmids
    tuple val(meta), path("*_plasmid_summary.filtered.tsv"), emit: summary_plasmids_filtered
    tuple val(meta), path("*_plasmid_ids.txt"), emit: sequence_ids
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Join plasmid summary tables
    head -1 ${plasmid_tables[0]} > ${prefix}_plasmids.tsv

    mkdir plasmids && mv *_plasmid_summary.tsv plasmids/ && cat plasmids/*_plasmid_summary.tsv | \\
        grep -v seq_name >> ${prefix}_plasmids.tsv && mv ${prefix}_plasmids.tsv ${prefix}_plasmid_summary.tsv

    # Filter plasmids using FDR and Hallmark genes:
    awk -F'\\t' 'BEGIN{OFS="\\t"}
    NR==1 { print; next } \\
    ( \\
    (\$2>=1000 && \$7<0.05 && \$3 ~ /DTR/) || \\
    (\$2>=1000 && \$7<0.05 && \$3 !~ /DTR/ && (\$8+0) >= 1) \\
    )' ${prefix}_plasmid_summary.tsv > ${prefix}_plasmid_summary.filtered.tsv

    cut -f1 ${prefix}_plasmid_summary.filtered.tsv > ${prefix}_plasmid_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version | head -n1 | cut -d' ' -f4))
    END_VERSIONS
    """
}
