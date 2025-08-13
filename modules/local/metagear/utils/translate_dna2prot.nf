process TRANSLATE_DNA2PROT {
    tag "$meta.id"
    label 'process_medium'
    // conda "conda-forge::python=3.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/auashen/biopython:1.83' :
        'docker.io/auashen/biopython:1.83' }"


    input:
    tuple val(meta), path(input_fp)

    output:
    tuple val(meta), path("*.prot.faa"), emit: prot_fasta_output
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def base = input_fp.baseName
    def output_fp = "${base}.prot.faa"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo running translate_fasta $input_fp $output_fp
    # python /nfs/data/work/shen/github/metagear-pipeline-internal/bin/
    translate_fasta.py $input_fp $output_fp
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: 3.8
    END_VERSIONS
    """

}
