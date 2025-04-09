process HUMANN_FUNCTION {
    maxForks 4
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::humann"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann:3.9--py312hdfd78af_0' :
        'biocontainers/humann:3.9--py312hdfd78af_0' }"

    input:
        tuple val(meta), path(fasta), path(profile)
        path humann3_uniref90
        path humann3_necleo

    output:
        tuple val(meta), path("${meta.id}/*_qc_genefamilies_cpm.tsv"), emit: gene_family
        tuple val(meta), path("${meta.id}/*_qc_pathabundance_cpm.tsv"), emit: path_abundance
        tuple val(meta), path("${meta.id}/*_qc_pathcoverage.tsv"), emit: path_coverage
        path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    gunzip -c ${fasta} > ${prefix}_qc.fasta

    humann --input ${prefix}_qc.fasta \\
        --input-format fasta --taxonomic-profile ${profile} \\
        --protein-database ${humann3_uniref90} \\
        --nucleotide-database ${humann3_necleo} \\
        -o ${prefix} \\
        --search-mode uniref90 \\
        --threads $task.cpus \\
        $args

    gzip ${prefix}_qc.fasta
    rm -r ${prefix}/${prefix}_qc_humann_temp

    humann_renorm_table --input ${prefix}/${prefix}_qc_genefamilies.tsv --output ${prefix}/${prefix}_qc_genefamilies_cpm.tsv $args2
    humann_renorm_table --input ${prefix}/${prefix}_qc_pathabundance.tsv --output ${prefix}/${prefix}_qc_pathabundance_cpm.tsv $args2

    rm ${prefix}/${prefix}_qc_genefamilies.tsv
    rm ${prefix}/${prefix}_qc_pathabundance.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann3: \$(humann --version 2>&1 | grep 'v3' | cut -d' ' -f2)
    END_VERSIONS
    """
}


process HUMANN_MERGE_PROFILES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::humann"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann:3.9--py312hdfd78af_0' :
        'biocontainers/humann:3.9--py312hdfd78af_0' }"

    input:
        tuple val(meta), path(profiles)

    output:
        tuple val(meta), path("*_merged_profiles.tsv"), emit: merged_profiles
        path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    humann_join_tables --input ./ --output ${prefix}_merged_profiles.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann3: \$(huamnn3 --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """
}