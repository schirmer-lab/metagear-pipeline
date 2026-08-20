#!/usr/bin/env bash
# Contract test for bin/vclust_to_pairs.py, which adapts Vclust's cluster output
# to the interface every consumer downstream of CLUSTER_SEQUENCES expects.
#
# Vclust writes `<sequence>\t<cluster_id>` keyed by a numeric cluster. MMseqs2's
# easy-cluster writes `representative\tmember` plus a representative FASTA, and
# that is the shape FIND_REPRESENTATIVES, the abundance subworkflow and the
# published catalog files are built on. If the conversion drifts, the vclust
# branch silently produces a catalog with a different interface from the mmseqs2
# branch, so the invariants are asserted here rather than left to a review.
#
# Asserted:
#   1. every input sequence appears exactly once as a member
#   2. one representative per cluster, and it is a LONGEST member of that
#      cluster (what CheckV's aniclust does), with the name breaking ties so the
#      choice is reproducible
#   3. each representative appears as its own member, the MMseqs2 convention
#   4. the representative FASTA holds exactly the representatives, no more
#   5. a sequence Vclust did not assign is a hard error, not a silent partial
#      catalog
#
# Runs the real script. No vclust installation needed: the cluster assignment is
# supplied directly, which is the script's actual input.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONVERT="$SCRIPT_DIR/../../bin/vclust_to_pairs.py"
[[ -f "$CONVERT" ]] || { echo "missing: $CONVERT" >&2; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

# Three clusters. c1 has a clear longest member (seq_long, 60 bp); c2 has two
# members of EQUAL length so the name tie-break is exercised; c3 is a singleton.
cat > in.fna <<'FASTA'
>seq_short
ACGTACGTAC
>seq_long
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>seq_mid
ACGTACGTACGTACGTACGT
>tie_b
ACGTACGTACGTACGTACGTACGTACGTACGT
>tie_a
ACGTACGTACGTACGTACGTACGTACGTACGT
>solo
ACGTACGTACGTACGTACGTACGTACGT
FASTA

cat > clusters.tsv <<'TSV'
object	cluster
seq_short	c1
seq_long	c1
seq_mid	c1
tie_b	c2
tie_a	c2
solo	c3
TSV

python "$CONVERT" --fasta in.fna --clusters clusters.tsv \
    --out-pairs pairs.tsv --out-representatives reps.fa.gz >/dev/null

# 1. every input sequence appears exactly once as a member
n_in="$(grep -c '^>' in.fna)"
n_members="$(cut -f2 pairs.tsv | wc -l)"
n_distinct="$(cut -f2 pairs.tsv | sort -u | wc -l)"
[[ "$n_members" == "$n_in" ]] || { echo "FAIL: $n_members member rows for $n_in sequences" >&2; exit 1; }
[[ "$n_distinct" == "$n_in" ]] || { echo "FAIL: members are not unique" >&2; exit 1; }

# 2a. the longest member represents c1
rep_c1="$(awk -F'\t' '$2=="seq_long"{print $1}' pairs.tsv)"
[[ "$rep_c1" == "seq_long" ]] || { echo "FAIL: c1 representative is $rep_c1, expected the longest member seq_long" >&2; exit 1; }
for member in seq_short seq_mid; do
    got="$(awk -F'\t' -v m="$member" '$2==m{print $1}' pairs.tsv)"
    [[ "$got" == "seq_long" ]] || { echo "FAIL: $member maps to $got, expected seq_long" >&2; exit 1; }
done

# 2b. equal lengths tie-break on name, so the result does not depend on input order
rep_c2="$(awk -F'\t' '$2=="tie_b"{print $1}' pairs.tsv)"
[[ "$rep_c2" == "tie_a" ]] || { echo "FAIL: equal-length tie gave $rep_c2, expected the lexicographically first name tie_a" >&2; exit 1; }

# 3. each representative is its own member
n_reps="$(cut -f1 pairs.tsv | sort -u | wc -l)"
n_selfrep="$(awk -F'\t' '$1==$2' pairs.tsv | wc -l)"
[[ "$n_reps" == "3" ]] || { echo "FAIL: $n_reps representatives, expected 3" >&2; exit 1; }
[[ "$n_selfrep" == "$n_reps" ]] || { echo "FAIL: $n_selfrep of $n_reps representatives are their own member" >&2; exit 1; }

# 4. the representative FASTA holds exactly the representatives
n_fa="$(zcat reps.fa.gz | grep -c '^>')"
[[ "$n_fa" == "$n_reps" ]] || { echo "FAIL: FASTA has $n_fa records for $n_reps representatives" >&2; exit 1; }
zcat reps.fa.gz | grep '^>' | sed 's/^>//' | sort > fa_names
cut -f1 pairs.tsv | sort -u > table_names
diff -q fa_names table_names >/dev/null || { echo "FAIL: FASTA names differ from table representatives" >&2; exit 1; }
# and the sequence travelled with the name
zcat reps.fa.gz | awk '/^>seq_long$/{getline; print length($0)}' > long_len
[[ "$(cat long_len)" == "60" ]] || { echo "FAIL: seq_long is $(cat long_len) bp in the FASTA, expected 60" >&2; exit 1; }

# 5. an unassigned sequence must fail loudly rather than yield a partial catalog
grep -v '^solo' clusters.tsv > partial.tsv
if python "$CONVERT" --fasta in.fna --clusters partial.tsv \
        --out-pairs bad.tsv --out-representatives bad.fa.gz >/dev/null 2>&1; then
    echo "FAIL: a missing sequence was accepted; a partial catalog would be published" >&2
    exit 1
fi

echo "PASS: test_vclust_to_pairs.sh"
