process PHABOX2_DOWNLOAD {
    tag "phabox"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/phabox:2.1.13--pyhdfd78af_1':
        'quay.io/biocontainers/phabox:2.1.13--pyhdfd78af_1' }"

    // Downloads the PhaBOX2 model and database bundle (~673 MB compressed).
    //
    // Fetched with wget rather than through the tool because phabox2 has NO
    // download subcommand: it is a single flat argument parser, and `phabox2
    // download` exits with "unrecognized arguments: download". The URL is
    // therefore version-bound and hardcoded, since there is no command to ask for
    // it. v2_2 is the current bundle and upstream directs users to it.
    //
    // The archive unpacks into a phabox_db_v2_2/ subdirectory alongside a
    // __MACOSX/ artefact, so the models are lifted up to a stable `phabox_db/`
    // path here. Consumers then get the same directory shape regardless of how
    // upstream names a future release.

    output:
    path "phabox_db/"   , emit: phabox_db
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def url = task.ext.url ?: 'https://github.com/KennthShang/PhaBOX/releases/download/v2/phabox_db_v2_2.zip'
    """
    wget --tries=5 --continue -O phabox_db.zip "${url}" ${args}
    unzip -q phabox_db.zip -d unpacked
    rm -f phabox_db.zip
    rm -rf unpacked/__MACOSX

    # Locate the directory that actually holds the models rather than assuming the
    # release's directory name, then promote it to a stable path. `bert/` is the
    # PhaTYP model directory and is the marker used to identify it.
    if [[ -d unpacked/bert ]]; then
        mv unpacked phabox_db
    else
        MODELDIR=""
        for d in unpacked/*/; do
            if [[ -d "\$d/bert" ]]; then MODELDIR="\$d"; break; fi
        done
        if [[ -z "\$MODELDIR" ]]; then
            echo "ERROR: no PhaBOX2 model directory (containing bert/) found in the archive." >&2
            echo "The release layout has changed; inspect the bundle at ${url}" >&2
            exit 1
        fi
        mv "\$MODELDIR" phabox_db
        rm -rf unpacked
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phabox: \$(phabox2 --help 2>&1 | sed -e 's/\\x1b\\[[0-9;]*m//g' | grep -om1 'PhaBOX v[0-9.]*' | sed 's/^PhaBOX v//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p phabox_db/bert
    touch phabox_db/bert/config.json
    touch phabox_db/bert/pytorch_model.bin
    mkdir -p phabox_db/database
    touch phabox_db/database/phatyp_db.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phabox: 2.1.13
    END_VERSIONS
    """
}
