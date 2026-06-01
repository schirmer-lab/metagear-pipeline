process MERGE_DREP_TABLES {
    tag "cohort"
    label 'process_single'

    // Pure-shell awk — same minimal container as STAGE_DREP_WORK.
    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5':
        'biocontainers/coreutils:9.5' }"

    // Concatenates per-batch dRep summary tables (Cdb.csv, Wdb.csv) into
    // cohort-level files. dRep emits these tables per invocation; since we
    // now run dRep K times (once per batch from BATCH_BINS), we need to
    // collapse them back into single tables for downstream analysis.
    //
    // Input files are pre-renamed by the upstream subworkflow to encode
    // the batch id in the filename:
    //   batch_<NNN>.Cdb.csv
    //   batch_<NNN>.Wdb.csv
    // We prepend a `batch_id` column derived from that prefix so cluster
    // IDs (which are batch-local in dRep) remain disambiguable across
    // batches in the merged cohort table.

    input:
    tuple val(meta), path(cdb_csvs, stageAs: 'cdb/*')
    path(wdb_csvs, stageAs: 'wdb/*')

    output:
    tuple val(meta), path("Cdb.cohort.csv"), emit: cdb
    tuple val(meta), path("Wdb.cohort.csv"), emit: wdb
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    # Per-batch Cdb / Wdb use BATCH-LOCAL cluster IDs (e.g. "1_1") that
    # collide across batches. To produce a globally consistent cohort
    # table — and so downstream consumers (e.g. SUMMARIZE_MAG_CATALOG)
    # can key on cluster IDs without manually splicing in batch_id — we
    # do two transformations per row:
    #
    #   1. Prepend a `batch_id` column (preserves provenance).
    #   2. Rewrite the cluster ID column itself so it's globally unique:
    #        Cdb.secondary_cluster + Cdb.primary_cluster → "<batch_id>.<id>"
    #        Wdb.cluster → "<batch_id>.<id>"
    #
    # The header is taken from the first input verbatim (after the
    # batch_id prefix) — the columns are unchanged, only the values in
    # the cluster columns are rewritten.
    #
    # Inputs are pre-renamed by the upstream subworkflow to
    # batch_<NNN>.<table>.csv so awk can recover the batch_id from the
    # filename without needing a sidecar metadata file.

    # ── Cdb merge — rewrite secondary_cluster (col 2) and primary_cluster
    #               (col 6, when present).
    merge_cdb() {
        local indir=\$1
        local outfile=\$2
        local first=1
        for f in \$(ls \${indir}/*.csv | sort); do
            local stem=\$(basename "\$f" .csv)
            local bid="\${stem%%.*}"
            if [[ \$first -eq 1 ]]; then
                IFS= read -r header < "\$f"
                printf 'batch_id,%s\\n' "\$header" > "\$outfile"
                first=0
            fi
            tail -n +2 "\$f" | awk -F',' -v OFS=',' -v bid="\$bid" '{
                # Identify the cluster columns by header position from the
                # dRep schema. Cdb is:
                #   genome,secondary_cluster,threshold,cluster_method,
                #   comparison_algorithm,primary_cluster
                # so secondary_cluster=\$2 and primary_cluster=\$6.
                \$2 = bid "." \$2
                if (NF >= 6) { \$6 = bid "." \$6 }
                print bid, \$0
            }' >> "\$outfile"
        done
    }

    # ── Wdb merge — rewrite cluster (col 2). Wdb schema:
    #               genome,cluster,score
    merge_wdb() {
        local indir=\$1
        local outfile=\$2
        local first=1
        for f in \$(ls \${indir}/*.csv | sort); do
            local stem=\$(basename "\$f" .csv)
            local bid="\${stem%%.*}"
            if [[ \$first -eq 1 ]]; then
                IFS= read -r header < "\$f"
                printf 'batch_id,%s\\n' "\$header" > "\$outfile"
                first=0
            fi
            tail -n +2 "\$f" | awk -F',' -v OFS=',' -v bid="\$bid" '{
                \$2 = bid "." \$2
                print bid, \$0
            }' >> "\$outfile"
        done
    }

    merge_cdb cdb Cdb.cohort.csv
    merge_wdb wdb Wdb.cohort.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(awk --version | head -1 | sed 's/^.*[Aa]wk //')
    END_VERSIONS
    """

    stub:
    """
    printf 'batch_id,genome,secondary_cluster\\n' > Cdb.cohort.csv
    printf 'batch_id,genome,cluster\\n' > Wdb.cohort.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: 9.5
    END_VERSIONS
    """
}
