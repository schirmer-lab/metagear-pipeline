process GTDBTK_DOWNLOAD_DB {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dram:1.5.0--pyhdfd78af_0' :
        'quay.io/biocontainers/dram:1.5.0--pyhdfd78af_0' }"

    output:
    path("release226_db/"), emit: database
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    # 1) Get a CA bundle via certifi
    python -m ensurepip -U >/dev/null 2>&1 || true
    python -m pip install --user -q certifi
    CERT="\$(python -c 'import certifi; print(certifi.where())')"

    # 2) Point *everything* at it
    export SSL_CERT_FILE="\$CERT" CURL_CA_BUNDLE="\$CERT" REQUESTS_CA_BUNDLE="\$CERT" GIT_SSL_CAINFO="\$CERT"

    parts="aa ab ac ad ae af ag ah ai aj ak al am an" && \\
    base="https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/auxillary_files/gtdbtk_package/split_package" && \\
    mkdir -p db && \\
    for p in \$parts; do echo "gtdbtk_r226_data.tar.gz.part_\$p"; done | \\
    xargs -n1 -P5 -I{} curl -fL --retry 5 -C - -O \$base/{} && \\
    cat gtdbtk_r226_data.tar.gz.part_* > gtdbtk_r226_data.tar.gz && rm -rf gtdbtk_r226_data.tar.gz.part_* && \\
    tar -xzf gtdbtk_r226_data.tar.gz -C db && mv db/release226 release226_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version 2>/dev/null) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """

}
