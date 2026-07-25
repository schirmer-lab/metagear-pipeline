// Fork of nf-core's mmseqs/databases. The query-side processes (MMSEQS_EASY_CLUSTER,
// MMSEQS_EASYTAXONOMY, MMSEQS_EASYSEARCH) all pin mmseqs2=17.b804f, so the conda
// env here matches. `mmseqs databases` also needs an HTTP client; bare mmseqs2
// containers don't ship one, so for the singularity/apptainer container we
// reuse nf-core's pre-built mmseqs2_wget multi-tool image (which carries
// mmseqs2=18.8cc5c — slightly newer than our query side; format-compatible
// for the GTDB taxonomy DB and only used for the one-time download).
process MMSEQS_DATABASES {
    tag "${database}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/edfecaaca16ca7fb7b6428dce0ed9c737549b38146360c98fdabf74e6c4cac68/data' :
        'community.wave.seqera.io/library/mmseqs2_wget:aa683a2c5355899d' }"

    input:
    val database

    output:
    path "${prefix}/"   , emit: database
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'mmseqs_database'
    """
    mkdir -p ${prefix}/

    mmseqs databases \\
        ${database} \\
        ${prefix}/database \\
        tmp/ \\
        --threads ${task.cpus} \\
        ${args}

    rm -rf tmp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: 'mmseqs_database'
    """
    mkdir -p ${prefix}/
    touch ${prefix}/database
    touch ${prefix}/database.dbtype
    touch ${prefix}/database_h
    touch ${prefix}/database_h.dbtype
    touch ${prefix}/database_h.index
    touch ${prefix}/database.index
    touch ${prefix}/database.lookup
    touch ${prefix}/database_mapping
    touch ${prefix}/database.source
    touch ${prefix}/database_taxonomy
    touch ${prefix}/database.version

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: 17.b804f
    END_VERSIONS
    """
}
