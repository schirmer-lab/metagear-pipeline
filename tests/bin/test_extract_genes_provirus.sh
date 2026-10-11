#!/usr/bin/env bash
# Regression test for provirus gene extraction in bin/extract_genes.sh.
#
# Before the fix, the awk block that assigns genes to provirus contigs by
# coordinate containment read the contig name from $1 with the leading FASTA
# ">" still attached: gsub(/^>/, "", header) rewrites a copy of the record, not
# $1 itself. The boundaries[] lookup therefore compared ">F-041_k141_225228"
# against a key of "F-041_k141_225228", never matched, and EVERY provirus gene
# was dropped. Cohort-wide effect when this was found: 17,610 (mgx), 41 (vmx)
# and 248 (mtx) provirus contigs in the per-sample viral id lists produced zero
# extracted genes, biased toward long, high-completeness prophages.
#
# The test drives the real script and asserts that a gene inside the provirus
# interval is kept while a gene beyond its end is not. It exercises the id path
# only, so seqtk is not required.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXTRACT_GENES="$SCRIPT_DIR/../../bin/extract_genes.sh"
[[ -f "$EXTRACT_GENES" ]] || { echo "missing: $EXTRACT_GENES" >&2; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

# One provirus representative: base contig F-041_k141_225228, interval 3-74691.
printf 'F-041_k141_225228|provirus_3_74691\n' > sequence_ids.txt

# Two genes on that base contig: one inside the interval, one past its end.
cat > genes.fasta <<'FASTA'
>F-041_k141_225228::2::100::500::+
ATGAAAGGGTTTTAA
>F-041_k141_225228::3::80000::80500::+
ATGCCCGGGTTTTAA
FASTA

# seqtk is only needed for the FASTA outputs; the id assertions below are what
# this test cares about, so a missing seqtk must not mask a real regression.
if ! command -v seqtk >/dev/null 2>&1; then
    printf '#!/usr/bin/env bash\ncat /dev/null\n' > seqtk
    chmod +x seqtk
    PATH="$WORK:$PATH"
    export PATH
fi

PREFIX=test_extract bash "$EXTRACT_GENES" sequence_ids.txt genes.fasta >/dev/null

PARTIAL="test_extract.partial.txt"
[[ -f "$PARTIAL" ]] || { echo "FAIL: no partial ids file produced" >&2; exit 1; }

inside=$(grep -cx 'F-041_k141_225228::2::100::500::+' "$PARTIAL" || true)
outside=$(grep -cx 'F-041_k141_225228::3::80000::80500::+' "$PARTIAL" || true)

status=0
if [[ "$inside" -ne 1 ]]; then
    echo "FAIL: gene inside the provirus interval (100-500) was not extracted" >&2
    status=1
fi
if [[ "$outside" -ne 0 ]]; then
    echo "FAIL: gene beyond the provirus interval (80000-80500) was extracted" >&2
    status=1
fi

if [[ "$status" -eq 0 ]]; then
    echo "PASS: provirus genes assigned by coordinate containment"
fi
exit "$status"
