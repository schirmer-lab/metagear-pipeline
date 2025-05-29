process SPLIT_FASTA {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/auashen/pyfasta:0.5.2' :
        'docker.io/auashen/pyfasta:0.5.2' }"

    input:
    tuple val(meta), path(input_fp), val(n_seq)

    output:
    path("splited_prot90"), emit: prot_fasta_split
    path "versions.yml", emit: versions


    script:
    def file_name = input_fp.getName()
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // pyfasta split -n 6 [-k 5000 ] some.fasta
    """
    mkdir -p splited_prot90;
    cp $input_fp splited_prot90/
    pyfasta split -n $n_seq splited_prot90/$file_name
    rm splited_prot90/$file_name
    rm splited_prot90/*.flat splited_prot90/*.gdx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pyfasta: \$(pip show pyfasta | grep "Version" | sed 's/Version: //g')
    END_VERSIONS
    """
}
