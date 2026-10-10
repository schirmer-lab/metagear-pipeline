// bwa reads the index and never the FASTA, so an index built from a different catalog
// maps to the wrong names without failing. This compares the sequence count in the
// index's .ann header against the catalog and passes the index through unchanged, so
// nothing downstream can use an index that was not built from this catalog.
process BWA_INDEX_CHECK {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/schirmer-lab/singularity-images/releases/download/23.11.22/python_3.10.sif' :
        'docker.io/raphsoft/python_base:3.10-R4' }"

    input:
    tuple val(meta), path(index), path(catalog)

    output:
    tuple val(meta), path(index), emit: index

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ann=\$(find -L ${index} -name '*.ann' | head -1)
    if [ -z "\$ann" ]; then
        echo "no .ann file under ${index}: not a bwa index" >&2
        exit 1
    fi

    read -r index_bases index_sequences _ < "\$ann"
    catalog_sequences=\$(grep -c '^>' ${catalog})

    if [ "\$index_sequences" != "\$catalog_sequences" ]; then
        echo "index was not built from ${catalog}: the index holds \$index_sequences" >&2
        echo "sequences over \$index_bases bases, the catalog holds \$catalog_sequences" >&2
        exit 1
    fi

    echo "index matches ${catalog}: \$index_sequences sequences over \$index_bases bases"
    """
}
