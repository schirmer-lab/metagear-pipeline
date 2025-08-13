process FUNCTIONALGROUP_ANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.27/python_3.10.sif' :
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    tuple val(meta), path(input_fp_lst)

    output:
    path("*.FG_IPS_Pfam.tsv"), emit: fg_ann_fp
    path("*.annotations.tsv"), emit: combined_annotations
    path "versions.yml", emit: versions

    script:
    def file_name = input_fp_lst.first().getName()
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // cat $input_fp_lst > interproscan_annotation_combined.tsv replaced into a better way
    """
    mkdir tmp_batches
    files=($input_fp_lst)

    # merge every 300 files to tmp_batches folder, then combine
    for ((i=0; i<\${#files[@]}; i+=300)); do
        batch_files="\${files[@]:i:300}"
        tmp_fp=tmp_batches/batch_\${i}.tsv
        cat \$batch_files > \$tmp_fp
    done
    cat tmp_batches/batch_*.tsv > ${prefix}.annotations.tsv

    functional_group_annotation.py ${prefix}.annotations.tsv ${prefix}.FG_IPS_Pfam.tsv
    rm -rf tmp_batches # remove the temp folder after running

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
