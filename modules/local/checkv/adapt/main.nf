process CHECKV_ADAPT_OUTPUT {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkv:1.0.1--pyhdfd78af_0':
        'biocontainers/checkv:1.0.1--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(original_file)

    output:
        tuple val(meta), path("${prefix}/*.${suffix}_*"), emit: adapted_results
        path "versions.yml"                             , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = task.ext.suffix ?: "${meta.suffix}"

    """
    mkdir -p ${prefix}

    # Resolve the ultimate target of the existing symlink
    target=\$(readlink -f "${original_file}")
    orig_filename=\$(basename "\$target")
    symlink_name="${prefix}/${prefix}.${suffix}_\${orig_filename}"

    # Create a new symlink pointing directly at the resolved file
    ln -s "\$target" "\${symlink_name}"

    printf '"%s":\n    checkv: %s\n' "${task.process}" "1.0.1" > versions.yml
    """
}
