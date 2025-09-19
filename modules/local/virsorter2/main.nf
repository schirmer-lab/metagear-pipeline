process VIRSORTER2 {
    tag "$meta.id"
    conda "bioconda::virsorter==2.2.4--pyhdfd78af_1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/virsorter:2.2.4--pyhdfd78af_0':
        'biocontainers/virsorter:2.2.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(contigs)
    path db

    output:
    tuple val(meta), path("*/*.final-viral-combined.fa"), emit: vs2_virus
    tuple val(meta), path("*/*.final-viral-score.tsv"), emit: vs2_score
    tuple val(meta), path("*/for-dramv/*.final-viral-combined-for-dramv.fa"), optional: true, emit: vs2_4dram_virus
    tuple val(meta), path("*/for-dramv/*.viral-affi-contigs-for-dramv.tab"), optional: true, emit: vs2_4dra_affi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: '' // Name for the resulting contig
    def vs2_folder = "${prefix}${args2}"
    """
    mkdir -p "\$PWD/.conda_pkgs"
    mkdir -p "\$PWD/_tmp"

    export CONDA_PKGS_DIRS="\$PWD/.conda_pkgs"
    export TMPDIR="\$PWD/_tmp" 

    INPUT=$contigs
    if [[ $contigs == *.gz ]]
    then
        gunzip -c $contigs > ${prefix}.fasta
        INPUT=${prefix}.fasta
    fi

    virsorter run $args --conda-prefix ./.conda -i \$INPUT -w ${vs2_folder} -j $task.cpus -d ${db} all

    mv ${vs2_folder}/final-viral-combined.fa ${vs2_folder}/${prefix}.final-viral-combined.fa
    mv ${vs2_folder}/final-viral-score.tsv ${vs2_folder}/${prefix}.final-viral-score.tsv
    mv ${vs2_folder}/for-dramv/final-viral-combined-for-dramv.fa ${vs2_folder}/for-dramv/${prefix}.final-viral-combined-for-dramv.fa || true
    mv ${vs2_folder}/for-dramv/viral-affi-contigs-for-dramv.tab ${vs2_folder}/for-dramv/${prefix}.viral-affi-contigs-for-dramv.tab || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        VirSorter2: \$(grep 'VirSorter' .command.log | head -n1 | sed 's/.* //')
    END_VERSIONS
    """
}
