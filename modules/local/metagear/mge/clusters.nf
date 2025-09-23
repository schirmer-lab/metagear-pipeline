process FIND_REPRESENTATIVES {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqtk=1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/schirmerlab/python310:25.09.10' :
        'docker.io/schirmerlab/python310:25.09.10' }"

    input:
    tuple val(meta), path(input_ids), path(input_clusters_tsv), path(input_representatives)

    output:
    tuple val(meta), path("*.representative_ids.txt"), emit: representative_ids
    tuple val(meta), path("*.representative.fa.gz"), emit: representatives
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -euo pipefail
    export LC_ALL=C

    IDS="${input_ids}"
    CLUSTERS="${input_clusters_tsv}"
    REPS_FA="${input_representatives}"

    PREFIX="${prefix}"
    OUT_IDS="\${PREFIX}.representative_ids.txt"
    OUT_FA="\${PREFIX}.representative.fa.gz"
    MISS_LOG="\${PREFIX}.missing_ids.txt"
    TMP_ALL="\${PREFIX}.rep_ids_in_order.tmp"
    TMP_IDS="\${PREFIX}.ids.txt"

    # Normalize (small) IDs file to plain text
    case "\$IDS" in
      *.gz) gzip -cd "\$IDS" > "\$TMP_IDS" ;;
      *)    cp -f  "\$IDS"    "\$TMP_IDS" ;;
    esac

    # Choose fastest available awk
    AWK_BIN="\$(command -v mawk || command -v gawk || command -v awk)"

    # Stream the (possibly huge) clusters file; do not materialize in temp
    if [[ "\$CLUSTERS" == *.gz ]]; then
      SRC_CMD=(gzip -cd "\$CLUSTERS")
    else
      SRC_CMD=(cat "\$CLUSTERS")
    fi

    # Map each query ID to its representative (first column).
    # - Loads only query IDs into RAM (O(|input_ids|)).
    # - Early-exits when all are found.
    # - Accepts tab or space as delimiter.
    # - If a queried ID is itself a representative, map to itself when seen in col1.
    "\${SRC_CMD[@]}" | "\$AWK_BIN" -v ids="\$TMP_IDS" -v miss="\$MISS_LOG" '
      BEGIN{
        FS="[[:space:]]+"; OFS="\\t";
        n=0; remaining=0;
        while ((getline line < ids) > 0) {
          if(!(line in want)) { want[line]=1; order[++n]=line; remaining++ }
        }
        close(ids)
      }
      {
        rep=\$1; mem=\$2
        if ((mem in want) && !(mem in found)) {
          found[mem]=rep; remaining--
        }
        if ((rep in want) && !(rep in found)) {
          found[rep]=rep; remaining--
        }
        if (remaining<=0) exit
      }
      END{
        for (i=1;i<=n;i++) {
          id=order[i]
          if (id in found) print found[id]
          else print id >> miss
        }
      }
    ' > "\$TMP_ALL"

    # Stable de-dup of representatives (preserve first-seen order)
    "\$AWK_BIN" '!seen[\$0]++' "\$TMP_ALL" > "\$OUT_IDS"

    # Extract sequences (seqtk supports gz input); ensure gzipped output
    if [[ -s "\$OUT_IDS" ]]; then
      seqtk subseq "\$REPS_FA" "\$OUT_IDS" | gzip -c > "\$OUT_FA"
    else
      : > "\${PREFIX}.empty.fa"
      gzip -c "\${PREFIX}.empty.fa" > "\$OUT_FA"
      rm -f "\${PREFIX}.empty.fa"
    fi

    # Optional summary (kept local to task workdir)
    {
      total=\$(wc -l < "\$TMP_IDS" | tr -d " ")
      mapped=\$(wc -l < "\$TMP_ALL" | tr -d " ")
      uniq=\$(wc -l < "\$OUT_IDS" | tr -d " ")
      echo "Queried: \$total"
      echo "Mapped (incl. dups): \$mapped"
      echo "Unique representatives: \$uniq"
      echo "Missing IDs listed in: \$MISS_LOG"
    } > "\${PREFIX}.summary.txt"

    # Cleanup small temps
    rm -f "\$TMP_ALL" "\$TMP_IDS"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk 2>&1 | sed -n 's/^Version: //p')
    END_VERSIONS
    """
}
