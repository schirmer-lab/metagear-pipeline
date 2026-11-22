#!/usr/bin/env bash
# Contract test for bin/msp_metaphlan_LM.py, which names each MSP after the
# MetaPhlAn species whose abundance fits it best across the cohort.
#
# The fit is only meaningful if an MSP abundance and a MetaPhlAn abundance are
# compared in the SAME sample. MSP_ABUNDANCE writes a header that labels the id
# column (`msp_name`), so the header has one field more than the values of a
# row. Reading that header as sample names pairs every value with the
# neighbouring sample and drops the last one: the script still runs, still
# writes R2 values, and names the wrong species.
#
# The fixture makes that failure visible. Each MSP tracks one species exactly,
# and a decoy species is built to fit perfectly under a one-sample shift.
#
# Asserted:
#   1. each MSP is matched to the species it tracks, with R2 of 1
#   2. every sample is used, none dropped
#   3. the same result for a header without an id label, and for one with an
#      empty id label
#   4. a row whose width does not fit the header is a hard error
#
# Runs the real script. Needs python with click, numpy and scikit-learn.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LM="$SCRIPT_DIR/../../bin/msp_metaphlan_LM.py"
[[ -f "$LM" ]] || { echo "missing: $LM" >&2; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

# msp_001 is 2 x species A, msp_002 is 3 x species B + 1.
printf 'msp_name\tS1\tS2\tS3\tS4\tS5\tS6\n'  > msp_labelled.tsv
printf 'msp_001\t2\t10\t4\t16\t6\t26\n'     >> msp_labelled.tsv
printf 'msp_002\t22\t4\t13\t1\t28\t10\n'    >> msp_labelled.tsv

# Sample columns in a different order from the MSP table, with the suffix the
# pipeline's merged profile carries. Decoy_A and Decoy_B hold the abundance of
# A and B in the NEXT sample, which is what a one-sample shift would compare.
{
    printf '#mpa_vTEST\n'
    printf 'clade_name\tS3_microbial\tS1_microbial\tS6_microbial\tS2_microbial\tS5_microbial\tS4_microbial\n'
    printf 'k__Bacteria\t100\t100\t100\t100\t100\t100\n'
    printf 'k__Bacteria|g__G|s__A\t2\t1\t13\t5\t3\t8\n'
    printf 'k__Bacteria|g__G|s__B\t4\t7\t3\t1\t9\t0\n'
    printf 'k__Bacteria|g__G|s__Decoy_A\t8\t5\t0\t2\t13\t3\n'
    printf 'k__Bacteria|g__G|s__Decoy_B\t0\t1\t0\t4\t3\t9\n'
} > metaphlan.tsv

best_species() { awk -F'\t' -v m="$2" '$1==m {n=split($2,a,"|"); print a[n]}' "$1"; }
best_r2()      { awk -F'\t' -v m="$2" '$1==m {print $6}' "$1"; }

check() {
    local label="$1" profile="$2" out="out_$1"
    python "$LM" --msp-profile "$profile" --metaphlan-profile metaphlan.tsv \
        --output-dir "$out" --metaphlan-version v4 > "$out.log"
    local best="$out/msp_metaphlan4_LM.bestR2.txt"
    # 1. the species each MSP tracks, at R2 of 1
    [[ "$(best_species "$best" msp_001)" == "s__A" ]] \
        || { echo "FAIL ($label): msp_001 matched $(best_species "$best" msp_001), expected s__A" >&2; exit 1; }
    [[ "$(best_species "$best" msp_002)" == "s__B" ]] \
        || { echo "FAIL ($label): msp_002 matched $(best_species "$best" msp_002), expected s__B" >&2; exit 1; }
    for m in msp_001 msp_002; do
        awk -v r="$(best_r2 "$best" "$m")" 'BEGIN { exit !(r > 0.999999) }' \
            || { echo "FAIL ($label): $m R2 is $(best_r2 "$best" "$m"), expected 1" >&2; exit 1; }
    done
    # 2. all six samples enter the fit
    grep -Eq '^N samples: +6$' "$out.log" \
        || { echo "FAIL ($label): $(grep '^N samples' "$out.log"), expected 6" >&2; exit 1; }
}

check labelled msp_labelled.tsv

# 3. the two other header conventions give the same answer
sed '1s/^msp_name\t//' msp_labelled.tsv > msp_unlabelled.tsv
sed '1s/^msp_name//'   msp_labelled.tsv > msp_empty_label.tsv
check unlabelled msp_unlabelled.tsv
check empty_label msp_empty_label.tsv

# 4. a row that does not fit the header must fail, not be zipped short
sed '2s/\t26$//' msp_labelled.tsv > msp_short_row.tsv
if python "$LM" --msp-profile msp_short_row.tsv --metaphlan-profile metaphlan.tsv \
        --output-dir out_bad --metaphlan-version v4 >/dev/null 2>&1; then
    echo "FAIL: a row narrower than the header was accepted" >&2
    exit 1
fi

echo "PASS: test_msp_metaphlan_lm.sh"
