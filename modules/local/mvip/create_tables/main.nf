process MERGE_TABLES {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.27/python_3.10.sif' :
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    tuple val(meta), val(files_with_meta), path(plasmid_summary)
    path ictv_taxonomy

    output:
    tuple val(meta), path("*_Merged_Genomad_CheckV_Summary.tsv"), emit: merged_tables, optional: true
    tuple val(meta), path("*_Merged_Genomad_CheckV_Summary.filtered.tsv"), emit: filtered_tables, optional: true
    tuple val(meta), path("*_viral_ids_to_keep.txt"), emit: sequence_ids, optional: true
    tuple val(meta), path("*_Plasmid_Summary.filtered.tsv"), emit: plasmid_filtered_tables, optional: true
    tuple val(meta), path("*_plasmid_ids_to_keep.txt"), emit: plasmid_sequence_ids, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def staged_files = files_with_meta.collect { file, pref ->
        def new_name = "${pref}.${file.getName()}"
        return "cp ${file} ${new_name}"
    }.join("\n")
    """
    $staged_files

    VIRUS_VIRUS_SUMMARY=\$(find . -maxdepth 1 -type f -name 'virus.${prefix}*.contigs_virus_summary.tsv' | head -n1)

    # Check if VIRUS_VIRUS_SUMMARY has more than 1 row (header + data)
    if [[ -f "\$VIRUS_VIRUS_SUMMARY" ]]; then
        VIRUS_ROWS=\$(wc -l < "\$VIRUS_VIRUS_SUMMARY")
        if [[ \$VIRUS_ROWS -le 1 ]]; then
            echo "No viral entries found in \$VIRUS_VIRUS_SUMMARY. Exiting."
            exit 0
        fi
    else
        echo "VIRUS_VIRUS_SUMMARY file not found. Exiting."
        exit 0
    fi

    merge_tables.py --sample-name ${prefix} \\
        --viral-checkv virus.quality_summary.tsv \\
        --provirus-checkv provirus.quality_summary.tsv \\
        --viral-genomad \$VIRUS_VIRUS_SUMMARY \\
        --provirus-genomad provirus.proviruses_virus_summary.tsv \\
        --ictv-taxonomy ${ictv_taxonomy} \\
        --output-file ./${prefix}_Merged_Genomad_CheckV_Summary.tsv

    # Create list of id to keep
    cat ./${prefix}_Merged_Genomad_CheckV_Summary.filtered.tsv | grep -v virus_id | cut -f2 > ${prefix}_viral_ids_to_keep.txt

    # Filter plasmids using FDR and Hallmark genes:
    awk -F'\\t' 'BEGIN{OFS="\\t"}
    NR==1 { print; next } \\
    ( \\
    (\$2>=1000 && \$7<0.05 && \$3 ~ /DTR/) || \\
    (\$2>=1000 && \$7<0.05 && \$3 !~ /DTR/ && (\$8+0) >= 1) \\
    )' ${plasmid_summary} > ${prefix}_Plasmid_Summary.filtered.tsv

    cut -f1 ${prefix}_Plasmid_Summary.filtered.tsv > ${prefix}_plasmid_ids_to_keep.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
