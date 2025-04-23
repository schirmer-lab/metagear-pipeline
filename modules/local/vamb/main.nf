process VAMB_CONCATENATE_FASTA {
    label 'process_medium'
    tag "${meta.id}"

    conda "bioconda::metaphlan=4.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://schirmerlab/tools/vamb:4.1.3' :
        'docker.io/schirmerlab/vamb:4.1.3' }"

    input:
    tuple val(meta), path(assemblies)

    output:
    tuple val(meta), path("*.fna.gz"), emit: catalog
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def args2 = task.ext.args2 ?: '.default' // Name for the resulting contig
    def m = args2 =~ /--catalog_name\s+(\S+)/
    def catalog_name = m.find() ? m.group(1) : ''

    """

    files=($assemblies)

    mkdir tmp_batches

    for ((i=0; i<\${#files[@]}; i+=1000)); do
        batch_files="\${files[@]:i:1000}"
        temp_file=./tmp_batches/temp_\${i}.fna.gz
        concatenate_fasta.py $args \$temp_file \$batch_files
        temp_files+=(\$temp_file)
    done

    concatenate_fasta.py $args ${prefix}.${catalog_name}.fna.gz \${temp_files[@]}

    # Remove the temporary files
    rm -rf tmp_batches

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vamb: \$(vamb --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}
