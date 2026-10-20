#!/usr/bin/env bash
# Runs the end-to-end test and checks its output. See docs/end_to_end_test.md.
set -euo pipefail

runtime=auto
outdir=""
keep=0

usage() {
    cat <<USAGE
Usage: tests/run_e2e.sh [--runtime singularity|apptainer|docker|conda|none] [--outdir DIR] [--keep]

Runs the gene-centric route to completion on the committed mock community and
checks the output. Takes a few minutes once the runtime has its images.

  --runtime  container runtime, or "none" to use tools already on PATH.
             Default: the first of singularity, apptainer, docker, conda found.
  --outdir   where to put results. Default: a temporary directory.
  --keep     keep the results and work directories on success.
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --runtime) runtime="${2:?--runtime needs a value}"; shift 2 ;;
        --outdir)  outdir="${2:?--outdir needs a value}";  shift 2 ;;
        --keep)    keep=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fail() { echo "FAIL: $*" >&2; exit 1; }

# Prerequisites
command -v java >/dev/null || fail "java not found. Nextflow needs a Java 17 or newer runtime, for example 'apt-get install openjdk-21-jre-headless'."
command -v nextflow >/dev/null || fail "nextflow not found. Install it with 'curl -s https://get.nextflow.io | bash' and put the resulting binary on PATH."

if [[ "$runtime" == auto ]]; then
    for candidate in singularity apptainer docker; do
        if command -v "$candidate" >/dev/null; then runtime="$candidate"; break; fi
    done
    if [[ "$runtime" == auto ]]; then
        for candidate in micromamba mamba conda; do
            if command -v "$candidate" >/dev/null; then runtime=conda; break; fi
        done
    fi
    [[ "$runtime" == auto ]] && fail "no runtime found. Install Apptainer or Singularity (preferred), or Docker, or a conda front end, then re-run. Pass --runtime none if the tools are already on PATH."
fi
case "$runtime" in
    docker) docker info >/dev/null 2>&1 || fail "the docker command is present but the daemon is not reachable. Start it, or add yourself to the docker group, or pass --runtime with another choice." ;;
    singularity|apptainer)
        "$runtime" --version >/dev/null 2>&1 || fail "$runtime is present but not working; check the installation."
        # Ubuntu 24.04 restricts unprivileged user namespaces, which a non-setuid
        # install needs. See docs/end_to_end_test.md.
        if [[ "$(cat /proc/sys/kernel/apparmor_restrict_unprivileged_userns 2>/dev/null || echo 0)" == 1 ]] \
           && ! compgen -G "/usr/lib*/apptainer/bin/starter-suid*" >/dev/null \
           && ! compgen -G "/usr/libexec/apptainer/bin/starter-suid*" >/dev/null; then
            echo "warning: unprivileged user namespaces are restricted on this kernel and no setuid starter was found."
            echo "         If the run fails with a user namespace error, install apptainer-suid. See docs/end_to_end_test.md."
        fi
        ;;
esac
echo "runtime: $runtime"

# Required for the pipeline to compile; harmless if already set.
export NXF_SYNTAX_PARSER=v1
# Where thread affinity is denied, as in some containers, MEGAHIT aborts without it.
export KMP_AFFINITY="${KMP_AFFINITY:-disabled}"

workroot="$(mktemp -d "${TMPDIR:-/tmp}/metagear_e2e.XXXXXX")"
[[ -n "$outdir" ]] || outdir="$workroot/results"
free_gb() { df -Pk "$1" | awk 'NR==2 {printf "%d", $4 / 1048576}'; }

if [[ "$runtime" == singularity || "$runtime" == apptainer ]]; then
    export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-$HOME/.singularity_cache}"
    mkdir -p "$NXF_SINGULARITY_CACHEDIR"
    echo "image cache: $NXF_SINGULARITY_CACHEDIR (first run downloads 6.6 GB, later runs reuse)"
    avail=$(free_gb "$NXF_SINGULARITY_CACHEDIR")
    [[ "$avail" -ge 15 ]] || echo "warning: only ${avail} GB free there; the images need about 15 GB. Set NXF_SINGULARITY_CACHEDIR somewhere larger."
fi
avail=$(free_gb "$workroot")
[[ "$avail" -ge 5 ]] || echo "warning: only ${avail} GB free in $workroot for the work directory. Set TMPDIR somewhere larger."

profile=test_e2e
[[ "$runtime" == none ]] || profile="test_e2e,$runtime"

started=$SECONDS
nextflow run "$repo" -profile "$profile" --outdir "$outdir" -work-dir "$workroot/work" || fail "the pipeline run did not complete; see the Nextflow error above"
elapsed=$((SECONDS - started))

# Checks. Counts vary between runs, so assert structure, never a sequence count.
checks=0
check() { checks=$((checks + 1)); echo "  ok: $1"; }

for f in mmseqs/all.genes.representative.fa.gz mmseqs/all.proteins.representative.fa.gz \
         mmseqs/all.genes.clusters.tsv coverm/all.genes_count_merged.tsv \
         coverm/all.genes_rpkm_merged.tsv coverm/all.genes_tpm_merged.tsv \
         pipeline_info/genes_mqc_versions.yml; do
    [[ -s "$outdir/$f" ]] || fail "missing or empty: $outdir/$f"
done
check "every expected output exists and is non-empty"

genes=$(zcat "$outdir/mmseqs/all.genes.representative.fa.gz" | grep -c '^>')
proteins=$(zcat "$outdir/mmseqs/all.proteins.representative.fa.gz" | grep -c '^>')
[[ "$genes" -gt 0 ]] || fail "the gene catalogue is empty"
[[ "$genes" -eq "$proteins" ]] || fail "gene catalogue has $genes sequences but the protein catalogue has $proteins"
check "gene and protein catalogues agree: $genes sequences each"

for matrix in count rpkm tpm; do
    f="$outdir/coverm/all.genes_${matrix}_merged.tsv"
    rows=$(($(wc -l < "$f") - 1))
    [[ "$rows" -eq "$genes" ]] || fail "$matrix matrix has $rows rows, expected $genes"
    cols=$(head -1 "$f" | tr '\t' '\n' | tail -n +2 | sort | tr '\n' ' ')
    [[ "$cols" == "MOCK1 MOCK2 MOCK3 " ]] || fail "$matrix matrix columns are '$cols', expected the three mock samples"
done
check "count, RPKM and TPM matrices are $genes genes by 3 samples"

for s in MOCK1 MOCK2 MOCK3; do
    [[ -s "$outdir/coverm/all.genes/$s.bam" ]] || fail "missing alignment: $outdir/coverm/all.genes/$s.bam"
done
check "all three samples aligned to the catalogue"

printf 'PASS: %s checks, %dm %02ds\n' "$checks" $((elapsed / 60)) $((elapsed % 60))
if [[ "$keep" -eq 1 ]]; then
    echo "results: $outdir"
    echo "work:    $workroot/work"
else
    rm -rf "$workroot"
fi
