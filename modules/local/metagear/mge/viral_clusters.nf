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
    tuple val(meta), path("*.representative.fa.gz"),  emit: representatives
    tuple val(meta), path("*.clusters.tsv"),          emit: representative_clusters
    tuple val(meta), path("*.pure_representatives.txt"), emit: pure_representatives
    tuple val(meta), path("*.clusters.annotated.tsv"), emit: input_clusters_annotated
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -euo pipefail
    export LC_ALL=C

    # Support multiple input_ids files (paths) by concatenating them first.
    printf '%s\\n' ${input_ids} | xargs cat >> all_input_ids.txt
    IDS=all_input_ids.txt

    CLUSTERS="${input_clusters_tsv}"
    REPS_FA="${input_representatives}"

    PREFIX="${prefix}"
    OUT_IDS="\${PREFIX}.representative_ids.txt"
    OUT_FA="\${PREFIX}.representative.fa.gz"
    OUT_CLU="\${PREFIX}.clusters.tsv"                 # cluster rows for selected representatives
    OUT_PURE="\${PREFIX}.pure_representatives.txt"
    OUT_ANN="\${PREFIX}.clusters.annotated.tsv"       # full clusters with annotation column
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

    # Stream the (possibly huge) clusters file; do not materialize
    if [[ "\$CLUSTERS" == *.gz ]]; then
      SRC_CMD=(gzip -cd "\$CLUSTERS")
    else
      SRC_CMD=(cat "\$CLUSTERS")
    fi

    #############################################
    # PASS 1: Map each input ID -> representative
    #############################################
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
        # member maps to its representative
        if ((mem in want) && !(mem in found)) {
          found[mem]=rep; remaining--
        }
        # if the representative itself is asked for, map to self
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

    #############################################
    # Extract representative sequences (gz in/out)
    #############################################
    if [[ -s "\$OUT_IDS" ]]; then
      seqtk subseq "\$REPS_FA" "\$OUT_IDS" | gzip -c > "\$OUT_FA"
    else
      : > "\${PREFIX}.empty.fa"
      gzip -c "\${PREFIX}.empty.fa" > "\$OUT_FA"
      rm -f "\${PREFIX}.empty.fa"
    fi

    #########################################################################
    # PASS 2: Single stream to:
    #  - write original cluster rows for selected representatives -> OUT_CLU
    #  - compute "pure" representatives (all members in input_ids)        -> OUT_PURE
    #  - write FULL clusters with an extra column = PREFIX iff member in input_ids -> OUT_ANN
    #########################################################################
    : > "\$OUT_PURE"
    : > "\$OUT_CLU"
    : > "\$OUT_ANN"

    "\${SRC_CMD[@]}" | "\$AWK_BIN" -v ids="\$TMP_IDS" -v repsf="\$OUT_IDS" -v out_tsv="\$OUT_CLU" -v out_pure="\$OUT_PURE" -v out_ann="\$OUT_ANN" -v label="\$PREFIX" '
      BEGIN{
        FS="[[:space:]]+"; OFS="\\t";
        # Set of input_ids
        while ((getline x < ids) > 0) want[x]=1
        close(ids)
        # Representatives of interest + order
        n=0
        while ((getline r < repsf) > 0) { reps[r]=1; reps_list[++n]=r }
        close(repsf)
        # Pre-mark reps not in want as mixed (in case no rep->rep line exists)
        for (i=1;i<=n;i++) { rep=reps_list[i]; if (!(rep in want)) mixed[rep]=1 }
      }
      {
        rep=\$1; mem=\$2

        # --- Full-annotation stream (for ALL rows) ---
        if (mem in want) print \$0, label >> out_ann
        else             print \$0        >> out_ann

        # --- Subset: rows of selected representatives + purity check ---
        if (rep in reps) {
          print \$0 >> out_tsv
          if (!(mem in want)) mixed[rep]=1
        }
      }
      END{
        close(out_ann)
        close(out_tsv)
        for (i=1; i<=n; i++) {
          rep=reps_list[i]
          if (!(rep in mixed)) print rep >> out_pure
        }
        close(out_pure)
      }
    '

    #############################################
    # Optional summary (for task logs)
    #############################################
    {
      total=\$(wc -l < "\$TMP_IDS" | tr -d " ")
      mapped=\$(wc -l < "\$TMP_ALL" | tr -d " ")
      uniq=\$(wc -l < "\$OUT_IDS" | tr -d " ")
      pure=\$(wc -l < "\$OUT_PURE" | tr -d " ")
      mixed=\$(( uniq - pure ))
      echo "Queried IDs:             \$total"
      echo "Mapped (incl. dups):     \$mapped"
      echo "Unique representatives:  \$uniq"
      echo "Pure clusters (all mems in input_ids): \$pure"
      echo "Mixed clusters:          \$mixed"
      echo "Missing IDs listed in:   \$MISS_LOG"
      echo "Annotated clusters TSV:  \$OUT_ANN"
    } > "\${PREFIX}.summary.txt"

    # Cleanup small temps
    rm -f "\$TMP_ALL" "\$TMP_IDS"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk 2>&1 | sed -n 's/^Version: //p')
    END_VERSIONS
    """
}


process MERGE_CLUSTER_ANNOTATIONS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqtk=1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/schirmerlab/python310:25.09.10' :
        'docker.io/schirmerlab/python310:25.09.10' }"

    input:
      tuple val(meta), path(annotated_clusters)

    output:
      tuple val(meta), path("all.genes.clusters.annotated.tsv"), emit: input_clusters_annotated
      path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -euo pipefail
    export LC_ALL=C

    # annotated_clusters may contain multiple files; merge column-3 labels by (col1,col2).
    # - Preserve first-seen row order across the concatenated inputs.
    # - Deduplicate labels per row, keep first-seen label order.
    # - Rows with no label in any file remain 2-column.
    AWK_BIN="\$(command -v mawk || command -v gawk || command -v awk)"

    # If any input is gz, transparently stream-decompress; otherwise cat.
    merge_stream() {
      for f in "\$@"; do
        case "\$f" in
          *.gz) gzip -cd "\$f" ;;
          *)    cat "\$f" ;;
        esac
      done
    }

    merge_stream ${annotated_clusters} | "\$AWK_BIN" '
      BEGIN{
        FS="[\\t ]+"; OFS="\\t";
        n=0;
      }
      {
        rep=\$1; mem=\$2;
        key = rep "\\t" mem;

        # First time we see this row, remember stable order
        if (!(key in seen_row)) {
          seen_row[key]=1;
          order[++n]=key;
        }

        # If a label exists in col3, aggregate it uniquely in first-seen order
        if (NF>=3 && \$3 != "") {
          lab=\$3;
          k = key SUBSEP lab;          # per-row-per-label uniqueness
          if (!(k in have)) {
            have[k]=1;
            if (labels[key]=="") labels[key]=lab;
            else labels[key]=labels[key] "," lab;
          }
        }
      }
      END{
        for (i=1; i<=n; i++) {
          key=order[i];
          if (key in labels) print key, labels[key];
          else                print key;
        }
      }
    ' > all.genes.clusters.annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk 2>&1 | sed -n 's/^Version: //p')
    END_VERSIONS
    """
}
