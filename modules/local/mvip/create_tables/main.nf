process CREATE_TABLES {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.27/python_3.10.sif' :
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    tuple val(meta), val(files_with_meta)
    path ictv_taxonomy

    output:
    tuple val(meta), path("*_Merged_Genomad_CheckV_Summary.tsv"), emit: merged_tables
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

    merge_tables.py --sample-name ${prefix} \\
        --virus-checkv virus.quality_summary.tsv \\
        --provirus-checkv provirus.quality_summary.tsv \\
        --virus-genomad virus.${prefix}.contigs_virus_summary.tsv \\
        --provirus-genomad provirus.proviruses_virus_summary.tsv \\
        --ictv-taxonomy ${ictv_taxonomy} \\
        --output-file ./${prefix}_Merged_Genomad_CheckV_Summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
