process MMSEQS_EASYTAXONOMY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:17.b804f--hd6d6fdc_1':
        'biocontainers/mmseqs2:17.b804f--hd6d6fdc_1' }"

    input:
    // contigs — per-sample input FASTA (e.g. the post-Binette unbinned set)
    // taxdb   — directory containing an mmseqs2 pre-built taxonomy DB (built
    //           via `mmseqs databases GTDB <out> <tmp>`); pass as a Path to
    //           the DB *prefix file* OR the containing directory (both
    //           supported by the resolver below).
    tuple val(meta), path(contigs)
    path taxdb

    output:
    tuple val(meta), path("${prefix}_lca.tsv")          , emit: lca           , optional: true
    tuple val(meta), path("${prefix}_tophit_report")    , emit: tophit_report , optional: true
    tuple val(meta), path("${prefix}_report")           , emit: report        , optional: true
    tuple val(meta), path("${prefix}")                  , emit: outdir
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    # Pre-flight memory check: mmseqs2 easy-taxonomy needs to fit (a chunk of)
    # the GTDB taxonomy DB into RAM. Recent GTDB builds in mmseqs2 format are
    # ~200 GiB on disk (was ~30 GiB historically, has grown ~6×). The
    # prefilter has ~14 GiB fixed overhead (query DB + working buffers); the
    # residual after that holds the k-mer index of the target-DB split chunk.
    # With the current GTDB build the smallest workable chunk needs >50 GiB,
    # so available memory must be at least ~96 GiB (split-memory-limit ~65 GiB
    # − 14 GiB overhead = 51 GiB residual).
    # /proc/meminfo's MemAvailable reflects any cgroup limit imposed by
    # params.max_memory + process.resourceLimits.
    AVAIL_MEM_G=\$(awk '/MemAvailable/ {printf "%d", \$2/1024/1024}' /proc/meminfo)
    if [ "\$AVAIL_MEM_G" -lt 96 ]; then
        echo "ERROR: MMSEQS_EASYTAXONOMY needs at least ~96 GiB of available memory" >&2
        echo "       to query the GTDB taxonomy DB. Currently seeing \${AVAIL_MEM_G} GiB." >&2
        echo "" >&2
        echo "Most likely cause: params.max_memory in your ~/.metagear/metagear.config is" >&2
        echo "set below the per-task memory this step needs." >&2
        echo "" >&2
        echo "Fix: raise params.max_memory to >=128.GB so this step gets its requested 128 GB." >&2
        echo "Nextflow's scheduler caps concurrency at (host_memory / task.memory) automatically;" >&2
        echo "no maxForks limit is needed." >&2
        exit 1
    fi

    mkdir -p ${prefix} tmp

    # `mmseqs databases GTDB <out> <tmp>` produces a DB whose handle is the
    # path prefix (e.g. /path/to/gtdb_database -> sibling files
    # gtdb_database.dbtype, .index, .lookup, gtdb_database_h, gtdb_database_h.dbtype,
    # gtdb_database_taxonomy, etc.).
    #
    # We accept either the DB prefix file directly OR a directory containing
    # the DB. For the directory case we use `.lookup` as the anchor — it's
    # unique to the main DB; companion DBs (database_h, database_taxonomy,
    # database_mapping, …) do NOT have a `.lookup` file, so this avoids
    # mistakenly handing the header DB to easy-taxonomy as the main handle
    # (which would surface as "Database <path>_h needs header information").
    if [ -d ${taxdb} ]; then
        LOOKUP=\$(find -L ${taxdb} -maxdepth 1 -name '*.lookup' -type f | head -n1)
        if [ -z "\$LOOKUP" ]; then
            echo "ERROR: ${taxdb} does not contain an mmseqs2 DB (no *.lookup found at depth 1)" >&2
            exit 1
        fi
        DB_PREFIX="\${LOOKUP%.lookup}"
    else
        DB_PREFIX=${taxdb}
    fi

    INPUT=${contigs}
    if [[ ${contigs} == *.gz ]]; then
        gunzip -c ${contigs} > input.fna
        INPUT=input.fna
    fi

    mmseqs easy-taxonomy \\
        \${INPUT} \\
        \${DB_PREFIX} \\
        ${prefix}/${prefix} \\
        tmp \\
        $args \\
        --threads ${task.cpus}

    # easy-taxonomy writes <prefix>_lca.tsv, <prefix>_tophit_report,
    # <prefix>_report inside the output dir; surface them at the top level for
    # easy emit globbing (keep the dir too for full evidence).
    for f in ${prefix}/${prefix}_lca.tsv ${prefix}/${prefix}_tophit_report ${prefix}/${prefix}_report; do
        if [ -f "\$f" ]; then
            cp "\$f" "\$(basename \$f)"
        fi
    done

    rm -rf tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}_lca.tsv ${prefix}_tophit_report ${prefix}_report
    touch ${prefix}/${prefix}_lca.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: 17.b804f
    END_VERSIONS
    """
}
