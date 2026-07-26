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
    tuple val(meta), path("*.faa.gz"), emit: prot_fasta_output
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // def base = input_fp.baseName
    // def output_fp = "${base}.faa"
    // def output_fp = base.endsWith(".fa") ? base.substring(0, base.length()-3) : "${base}.faa"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id.replace(".genes", ".proteins")}"

    """

    INPUT=$input_fp
    if [[ $input_fp == *.gz ]]
    then
        gunzip -c $input_fp >| \$PWD/${prefix}_plain.fa
        INPUT=\$PWD/${prefix}_plain.fa
    fi

    echo running translate_fasta \$INPUT "${prefix}.faa"
    translate_fasta.py \$INPUT "${prefix}.faa"

    gzip "${prefix}.faa"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: 3.8
    END_VERSIONS
    """

}
