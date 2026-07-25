process PHOLD_INSTALL {
    tag "phold"
    label 'process_high'

    conda "bioconda::phold=1.2.5"
    // quay.io only — Galaxy-depot's SIF mirror is unreliable for this image
    // (~1.1 GB; routinely drops mid-pull). Singularity pulls from quay.io via
    // docker:// and converts to SIF locally with no Docker runtime needed.
    container 'quay.io/biocontainers/phold:1.2.5--pyhdfd78af_0'

    // Downloads and indexes the PHOLD structural reference DB.
    //
    // Default DB (~7.7 GB on disk): ~1.36M Foldseek-indexed structures
    // derived from PHROG representatives. Extended DB (~9–10 GB) adds 1.8M
    // efam + enVhog proteins as "ghost" matches; opt in via
    //   ext.args = '--medium'
    // in the workflow-level config when desired.
    //
    // The `--foldseek_gpu` flag (also via ext.args) builds a GPU-accelerated
    // Foldseek index; only enable on clusters where compare-time runs on
    // GPU nodes.

    output:
        path "phold_db/"    , emit: phold_db
        path "versions.yml" , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Force huggingface_hub's legacy Python downloader. The default Rust-based
    # hf-xet downloader (huggingface_hub >= 0.30) does its own CA-bundle lookup
    # via reqwest, which panics with "No CA certificates were loaded from the
    # system" inside the PHOLD biocontainer (no /etc/ssl/certs shipped). On
    # panic, phold falls back to Zenodo — and v1.2.5's Zenodo path has a
    # signature bug (download_requests called with 4 args but takes 2). So this
    # opt-out env var keeps us on the HF happy-path and avoids the broken
    # fallback. Drop once phold's container ships CA roots or the fallback bug
    # is fixed upstream.
    export HF_HUB_DISABLE_XET=1

    phold install \\
        -d phold_db \\
        -t ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phold: \$(phold --version 2>&1 | sed 's/^phold, version //; s/^phold //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p phold_db
    touch phold_db/phold_structure_foldseek_db
    touch phold_db/phold_structure_foldseek_db.index
    touch phold_db/phold_annotations.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phold: 1.2.5
    END_VERSIONS
    """
}
