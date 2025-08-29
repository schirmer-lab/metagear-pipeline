process DRAM_SETUP {
    label 'process_single'

    conda "bioconda::dram==1.5.0--pyhdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dram:1.5.0--pyhdfd78af_0' :
        'quay.io/biocontainers/dram:1.5.0--pyhdfd78af_0' }"

    output:
    path "dram_db/", emit: dram_db
    path "versions.yml", emit: versions

    script:
    """
    mkdir dram_db
    export TMPDIR="\$PWD/.tmp"; mkdir -p "\$TMPDIR"
    export DRAM_CONFIG_LOCATION=./dram_db/CONFIG

    cat > dram_db/CONFIG <<'JSON'
    {
        "vogdb": "",
        "vog_annotations": "",
        "uniref": null
    }
JSON

    # 1) Get a CA bundle via certifi
    python -m ensurepip -U >/dev/null 2>&1 || true
    python -m pip install --user -q certifi
    CERT="\$(python -c 'import certifi; print(certifi.where())')"

    # 2) Point *everything* at it
    export SSL_CERT_FILE="\$CERT" CURL_CA_BUNDLE="\$CERT" REQUESTS_CA_BUNDLE="\$CERT" GIT_SSL_CAINFO="\$CERT"

    # 3) Tell wget explicitly (robust across OpenSSL/GnuTLS builds)
    printf 'ca_certificate = %s\ncheck_certificate = on\n' "\$CERT" > .wgetrc
    export WGETRC="\$PWD/.wgetrc"

    # 4) Download VOG separately to avoid "latest" pointer failures
    VOGREL=231
    mkdir -p vog_seed

    echo "Downloading VOG databases (version \${VOGREL} - final version)..."
    curl -fL -C - --retry 3 -o vog_seed/vog.hmm.tar.gz "https://fileshare.lisc.univie.ac.at//vog/vog\${VOGREL}/vog.hmm.tar.gz"
    curl -fL -C - --retry 3 -o vog_seed/vog.annotations.tsv.gz "https://fileshare.lisc.univie.ac.at//vog/vog\${VOGREL}/vog.annotations.tsv.gz"

    # 4.1) Re-package VOG database (remove hmm/ folder so DRAM-setup.py can use it)
    cwd1="\$PWD" && cd vog_seed && workdir=".vog_hmm_work_\$\$" && mkdir -p "\$workdir" && \\
    tar -xzf vog.hmm.tar.gz -C "\$workdir" && \\
    find "\$workdir" -mindepth 2 -type f -name '*.hmm' -exec mv -f {} "\$workdir"/ \\; && \\
    cwd="\$PWD" && cd "\$workdir" && tar -cf - -- *.hmm | gzip -1 > "\$cwd"/vog.hmm.flat.tar.gz && cd "\$cwd" && \\
    rm -rf "\$workdir" && cd "\$cwd1"

    # 5) Download Pfam files separately to avoid hanging
    mkdir -p pfam_seed
    echo "Downloading Pfam-A.full.gz..."
    curl -fL -C - --retry 3 --max-time 3600 -o pfam_seed/Pfam-A.full.gz "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz" || \\
    curl -fL -C - --retry 3 --max-time 3600 -o pfam_seed/Pfam-A.full.gz "http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz"

    echo "Downloading Pfam-A.hmm.dat.gz..."
    curl -fL -C - --retry 3 --max-time 1800 -o pfam_seed/Pfam-A.hmm.dat.gz "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz" || \\
    curl -fL -C - --retry 3 --max-time 1800 -o pfam_seed/Pfam-A.hmm.dat.gz "http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz"

    # 6) Run DRAM setup
    echo "Running DRAM setup"
    DRAM-setup.py prepare_databases \\
        --output_dir dram_db \\
        --threads ${task.cpus} \\
        --skip_uniref \\
        --pfam_loc pfam_seed/Pfam-A.full.gz \\
        --pfam_hmm_loc pfam_seed/Pfam-A.hmm.dat.gz \\
        --vogdb_loc vog_seed/vog.hmm.flat.tar.gz \\
        --vog_annotations vog_seed/vog.annotations.tsv.gz \\
        --verbose | tee -a dram_db/database_processing.log

    # 7) Update paths in DRAM config
    sed -E -i.bak 's#"/[^"]*/dram_db/#"${params.dram_db}/#g' dram_db/CONFIG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        DRAM: 1.5.0
    END_VERSIONS
    """
}
