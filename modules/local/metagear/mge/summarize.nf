process COLLECT_TABLES {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/25.08.11/python_3.10.sif':
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    tuple val(meta), path(merged_virus_tables), path(plasmid_tables)

    output:
    tuple val(meta), path("*_checkv_summary_taxa.tsv"), emit: summary_taxa
    tuple val(meta), path("*_plasmid_summary.tsv"), emit: summary_plasmids
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Join plasmid summary tables
    head -1 ${merged_virus_tables[0]} > ${prefix}_checkv_summary_taxa.tsv

    mkdir virus && mv *_Merged_Genomad_CheckV_Summary.tsv virus/ && cat virus/*_Merged_Genomad_CheckV_Summary.tsv | \\
        grep -v virus_id >> ${prefix}_checkv_summary_taxa.tsv

    # Join plasmid summary tables
    head -1 ${plasmid_tables[0]} > ${prefix}_plasmids.tsv

    mkdir plasmids && mv *_plasmid_summary.tsv plasmids/ && cat plasmids/*_plasmid_summary.tsv | \\
        grep -v seq_name >> ${prefix}_plasmids.tsv && mv ${prefix}_plasmids.tsv ${prefix}_plasmid_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version | head -n1 | cut -d' ' -f4))
    END_VERSIONS
    """
}
