process MERGE_PHOLD_PREDICTIONS {
    tag "cohort"
    label 'process_single'

    // python3 (rather than the minimal coreutils container we use for some
    // other merge steps) — the JSON file emitted per shard is a dict keyed by
    // protein_id, so a text-cat would produce invalid JSON. python's json
    // module is the cleanest way to dict-merge across shards.
    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11':
        'biocontainers/python:3.11' }"

    // PHOLD_PREDICT runs per-shard (SEQKIT_SPLIT2 chunks). Each shard emits
    // a predict_dir/ containing phold_aa.fasta + phold_3di.fasta. PHOLD's
    // proteins-compare expects ONE prediction_dir, not many — so we
    // concatenate the per-shard FASTAs into a single directory the
    // PHOLD_COMPARE process can consume.
    //
    // Concatenation is safe: every protein has a unique ID across shards
    // (the subset FASTA from FILTER_STRUCTURES_INPUTS is partitioned by
    // SEQKIT_SPLIT2, no duplicates).

    input:
    tuple val(meta), path(predict_dirs, stageAs: 'predict_dirs/*')

    output:
    tuple val(meta), path("merged_predict_dir/"), emit: merged_predict_dir
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail
    mkdir -p merged_predict_dir

    # PHOLD's proteins-predict emits FIVE files per shard. We merge the four
    # that proteins-compare downstream consumes (skipping the .log file):
    #
    #   phold_aa.fasta                           — AA sequences (cat)
    #   phold_3di.fasta                          — 3Di-encoded sequences (cat)
    #   phold_prostT5_3di_mean_probabilities.csv — per-protein confidence (cat, header-aware)
    #   phold_prostT5_3di_all_probabilities.json — per-residue probabilities (dict-merge in python)
    #
    # Concatenation is safe because every protein has a unique ID across
    # shards — the input FASTA is partitioned by SEQKIT_SPLIT2 with no
    # duplicates.

    # FASTAs — plain cat
    for f in predict_dirs/*/phold_aa.fasta; do
        [[ -f "\$f" ]] && cat "\$f" >> merged_predict_dir/phold_aa.fasta
    done

    for f in predict_dirs/*/phold_3di.fasta; do
        [[ -f "\$f" ]] && cat "\$f" >> merged_predict_dir/phold_3di.fasta
    done

    # CSV — keep header from first shard, skip on subsequent
    first=1
    for f in predict_dirs/*/phold_prostT5_3di_mean_probabilities.csv; do
        [[ -f "\$f" ]] || continue
        if [[ \$first -eq 1 ]]; then
            cat "\$f" > merged_predict_dir/phold_prostT5_3di_mean_probabilities.csv
            first=0
        else
            tail -n +2 "\$f" >> merged_predict_dir/phold_prostT5_3di_mean_probabilities.csv
        fi
    done

    # JSON — each shard's file is a dict keyed by protein_id; merge via python
    # since text-cat would produce invalid JSON. The heredoc body is left-
    # justified so python's parser doesn't see unexpected leading indentation.
python3 <<'PY'
import json, glob, sys
merged = {}
for path in sorted(glob.glob('predict_dirs/*/phold_prostT5_3di_all_probabilities.json')):
    with open(path) as fh:
        try:
            shard = json.load(fh)
        except json.JSONDecodeError as e:
            sys.stderr.write(f"WARN: skipping unreadable JSON {path}: {e}\\n")
            continue
    if isinstance(shard, dict):
        # Protein IDs are unique across shards (SEQKIT_SPLIT2 partitions the
        # input), so a plain update is safe; flag if any collide.
        collisions = set(shard) & set(merged)
        if collisions:
            sys.stderr.write(f"WARN: {len(collisions)} duplicate protein IDs in {path}\\n")
        merged.update(shard)
    else:
        sys.stderr.write(f"WARN: unexpected JSON shape (not a dict) in {path}\\n")
with open('merged_predict_dir/phold_prostT5_3di_all_probabilities.json', 'w') as out:
    json.dump(merged, out)
print(f"merged {len(merged)} protein-level probability entries", file=sys.stderr)
PY

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python3 --version 2>&1 | sed 's/Python //')
END_VERSIONS
    """

    stub:
    """
    mkdir -p merged_predict_dir
    touch merged_predict_dir/phold_aa.fasta
    touch merged_predict_dir/phold_3di.fasta
    touch merged_predict_dir/phold_prostT5_3di_mean_probabilities.csv
    echo '{}' > merged_predict_dir/phold_prostT5_3di_all_probabilities.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
