process PHOLD_PREDICT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::phold=1.2.5"
    // quay.io only — Galaxy-depot's SIF mirror is unreliable for this image
    // (~1.1 GB; routinely drops mid-pull). Singularity pulls from quay.io via
    // docker:// and converts to SIF locally with no Docker runtime needed.
    container 'quay.io/biocontainers/phold:1.2.5--pyhdfd78af_0'

    // ProstT5 inference: protein AA sequences → Foldseek's 3Di alphabet.
    // This is the heavy step (ProstT5 is a 3B-param T5 fine-tune). Runs on
    // a SHARD of the input FASTA — the structures subworkflow splits the
    // catalog via SEQKIT_SPLIT2 before fanning out across PHOLD_PREDICT.
    //
    // GPU vs CPU is controlled by ext.args from conf/metagear/structures.config
    // (set `--cpu` when running on a non-GPU node). ProstT5 batch size is
    // tunable via ext.args (`--batch_size N`); default 1.
    //
    // PHOLD's prediction step accepts a protein FASTA via `proteins-predict`
    // and emits a directory containing phold_aa.fasta and phold_3di.fasta.
    // We pass meta.id (the shard tag) into ext.prefix so downstream merge
    // can deduplicate / order shards predictably.

    input:
    tuple val(meta), path(proteins_fasta)
    path phold_db

    output:
    tuple val(meta), path("predict_${meta.id}/"), emit: predict_dir
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Disable hf-xet (huggingface_hub's Rust downloader) — it panics with
    # "No CA certificates were loaded from the system" inside the phold
    # biocontainer. proteins-predict loads ProstT5 from the local phold DB
    # so HF shouldn't be touched, but if the model is missing locally PHOLD
    # tries to re-download, and that path is broken without this opt-out.
    # See modules/local/phold/install/main.nf for the full story.
    export HF_HUB_DISABLE_XET=1

    # `-d \${phold_db}` is required even for proteins-predict: PHOLD validates
    # the DB on startup (looks for acrs_plddt_over_70_metadata.tsv etc.) even
    # though the actual prediction only uses ProstT5 weights from the same dir.
    phold proteins-predict \\
        -i ${proteins_fasta} \\
        -o predict_${meta.id} \\
        -d ${phold_db} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phold: \$(phold --version 2>&1 | sed 's/^phold, version //; s/^phold //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p predict_${meta.id}
    touch predict_${meta.id}/phold_aa.fasta
    touch predict_${meta.id}/phold_3di.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phold: 1.2.5
    END_VERSIONS
    """
}
