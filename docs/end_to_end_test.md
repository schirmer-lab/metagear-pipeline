# End-to-end test

The `test` profile checks that the pipeline installs, that its configuration loads and that containers resolve. It dispatches to no analysis workflow, because `params.workflow = "test"` matches none of the branches in `workflows/metagear.nf`, so a green `test` run says nothing about whether an analysis produces output.

The `test_e2e` profile fills that gap. It runs the gene-centric route to completion on a committed mock community, with no reference database and, once the container images are present, no network access.

## Before anything: the supported way to run MetaGEAR is the command line interface

For real work, install and run through `metagear-tools`, which pins a pipeline release, writes a host configuration and assembles the per-workflow configuration at invocation:

```bash
git clone https://github.com/schirmer-lab/metagear-tools && cd metagear-tools && ./install.sh
metagear genes --input samples.csv
```

The profiles below exist to test the pipeline itself, in continuous integration and during development, where no wrapper is present. They are not the recommended way to analyse data.

## Running it

```bash
tests/run_e2e.sh
```

That is the whole thing. The script picks a runtime, sets the two environment variables the run needs, executes the profile and then checks the output, printing one line per check and the wall clock. It exits non-zero if anything is missing, so it is safe to run unattended.

Useful options:

```bash
tests/run_e2e.sh --runtime singularity   # or apptainer, docker, conda, or none for tools already on PATH
tests/run_e2e.sh --outdir ./e2e --keep   # keep the results and work directories to look at
```

Without `--runtime` the script takes the first of Singularity, Apptainer or Docker that it finds, and falls back to the conda profile if a conda front end is present. Results go to a temporary directory and are removed on success unless `--keep` is given.

The equivalent by hand, if you would rather see the Nextflow invocation:

```bash
export NXF_SYNTAX_PARSER=v1
nextflow run . -profile test_e2e,singularity --outdir results
```

Via nf-test, where the `e2e` tag selects it on its own:

```bash
nf-test test --tag e2e --profile +singularity
```

The test file sets `profile "test_e2e"`, which overrides `profile = "test"` in `nf-test.config`; without it the suite runs the plumbing profile and the assertions fail. The `+` form appends, so the command above resolves to `-profile test_e2e,singularity`. Continuous integration passes no `--tag`, so this test also runs on every pull request, once per profile in the `nf-test.yml` matrix.

Two environment variables matter, and the script exports both. `NXF_SYNTAX_PARSER=v1` is required for the pipeline to compile at all, not just here. `KMP_AFFINITY=disabled` is needed wherever setting thread affinity is denied, as in some container and sandbox configurations, because MEGAHIT otherwise aborts with `OMP: Error #179: Function pthread_setaffinity_np() failed` and exit code 250 on any input.

## From a clean machine

Written for Ubuntu 24.04 LTS (noble), which is what we test on. Nothing here is pipeline-specific: a Java runtime, Nextflow, a container runtime and `git`. As root or with `sudo`:

```bash
# 1. Base tools and Java. A minimal cloud image has neither curl nor git.
sudo apt-get update
sudo apt-get install -y curl git software-properties-common openjdk-21-jre-headless

# 2. Apptainer, from its official PPA. The -suid package is the one to use on
#    noble, which restricts unprivileged user namespaces by default.
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt-get update
sudo apt-get install -y apptainer apptainer-suid

# 3. Nextflow.
curl -s https://get.nextflow.io | bash
sudo install -m 0755 nextflow /usr/local/bin/nextflow

# 4. The pipeline, and the test.
git clone https://github.com/schirmer-lab/metagear-pipeline
cd metagear-pipeline
tests/run_e2e.sh
```

Confirm the prerequisites before step 4 if you want to fail early:

```bash
java -version && nextflow -v && apptainer --version
```

`openjdk-21-jre-headless` comes from noble's own repositories, so Java needs no PPA. Nextflow 26.x requires Java 17 or newer and 21 is what noble gives you by default. For other distributions, follow the Apptainer admin guide: the PPA is Ubuntu only, Debian uses the release `.deb` and Red Hat derivatives take it from EPEL.

### The noble trap worth knowing before you start

Ubuntu 24.04 sets `kernel.apparmor_restrict_unprivileged_userns=1`, which blocks the unprivileged user namespaces that a non-setuid Apptainer needs. Installing `apptainer-suid` alongside `apptainer`, as in step 2 above, avoids the problem by not needing them. If you install only `apptainer` and it fails with a user namespace or mount namespace error, that restriction is the cause; `sudo sysctl -w kernel.apparmor_restrict_unprivileged_userns=0` confirms it, and the Apptainer admin guide covers the durable options.

### If you would rather use Docker

```bash
curl -fsSL https://get.docker.com | sudo sh
sudo usermod -aG docker "$USER"    # then log out and back in, or use sudo
tests/run_e2e.sh --runtime docker
```

The script refuses to start if the Docker daemon is unreachable, which is what you get before the group change takes effect.

### Without root

If you cannot install system packages, install Java and Nextflow into a conda environment with micromamba and run `tests/run_e2e.sh --runtime conda`, which makes Nextflow build a small environment per process instead of using containers. That is slower on the first run and is the path continuous integration exercises for the conda profile. Apptainer itself needs either a package install or an unprivileged binary installation, so there is no container route that needs nothing from the system.

### Disk

The image cache needs about 15 GB. The script puts it in `$HOME/.singularity_cache` unless `NXF_SINGULARITY_CACHEDIR` is already set, warns if that filesystem has less than 15 GB free, and keeps the Nextflow work directory under `$TMPDIR`, which on a small root volume is worth pointing somewhere larger.

### How long it takes

The computation is about two minutes. Four timed runs of the analysis itself, with the tools already present, took 2 min 09 s, 2 min 07 s, 1 min 51 s and 1 min 32 s on 8 CPUs, and a full `tests/run_e2e.sh` invocation including its checks took 1 min 42 s.

Everything else is one-off setup, and the container images dominate it. The route pulls 11 distinct images, 6.62 GB in total as measured from the registries:

| image | used by | size |
| --- | --- | --- |
| `schirmerlab/vamb:4.1.3` | catalogue concatenation | 4,313 MB |
| `python_3.10.sif` (23.11.27) | gene filtering | 527 MB |
| `python_3.10.sif` (23.11.22) | samplesheet check, file renaming | 523 MB |
| MultiQC (Seqera) | summary | 459 MB |
| `schirmerlab/coverm_bwamem2:0.7.0` | read mapping | 194 MB |
| `mmseqs2:17.b804f` | clustering | 164 MB |
| `auashen/biopython:1.83` | translation | 137 MB |
| `coverm:0.7.0` | abundance | 134 MB |
| BWA (Seqera) | indexing | 92 MB |
| `megahit:1.2.9` | assembly | 72 MB |
| Prodigal (mulled) | gene calling | 7 MB |

Measured on a clean Ubuntu 24.04 VM following the steps above, with an empty image cache: **8 min 31 s** for the pipeline and 8 min 40 s including the script's checks, all 29 tasks successful. That is the number to expect for a first run. The second run on the same machine is just the compute, since the images are cached.

For reference, 6.62 GB at link speed is 9 minutes on 100 Mbit/s, 2 minutes on 500 Mbit/s and 1 minute on 1 Gbit/s, and Apptainer also has to convert the three Docker Hub images into SIFs locally, the 4.31 GB VAMB image among them.

One image accounts for 65% of that download. `VAMB_CONCATENATE_FASTA` pulls 4.3 GB of VAMB for a single call to `vamb.vambtools.concatenate_fasta`, which concatenates FASTA files and renames their headers. Replacing that call would remove the largest single cost in both this test and a first production run.

## What it covers

Input validation, assembly (MEGAHIT), gene calling (Prodigal), the completeness filter, gene clustering (MMseqs2), protein calling and protein clustering, read mapping (BWA) and abundance (CoverM) in counts, RPKM and TPM, plus the version report and MultiQC summary.

A successful run produces a gene catalogue of a few dozen sequences, a protein catalogue of the same size, and an abundance matrix with one row per representative and one column per mock sample.

The exact count varies between runs on identical input: nine runs gave between 64 and 75 representatives (64, 65, 66, 66, 67, 69, 72, 72, 75), from 124 to 133 called genes. MEGAHIT returns the same contigs each time, with the same lengths, but numbers and orders them differently, and that order propagates through gene calling and cluster representative selection. The assertions therefore check structure, that each stage produced a non-empty output with the expected sample columns, and never a sequence count. Do not add a count assertion here.

## The warnings a direct run prints

A direct `nextflow run` of a workflow that reads a parameter the wrapper normally supplies prints `Access to undefined parameter` for it. That is harmless: an undeclared parameter reads as null, null is falsey, and the branch is skipped.

The eight reuse parameters (`contigs_dir`, `genes_dir`, `representative_genes`, `representative_proteins`, `representative_proteins_annotations`, `representative_genes_tpm`, `representative_genes_rpkm`, `representative_genes_count`) are declared in `nextflow.config` and `nextflow_schema.json` with the same `false` defaults that `conf/metagear/common.config` sets, so neither a direct run nor a wrapper run warns about them. The remaining per-workflow inputs and database parameters are still undeclared by design, so a direct run of `virus`, `classification`, `mag` or `msp` will warn about whichever of those it reads.

A second warning, `The operator 'first' is useless when applied to a value channel`, comes from the `.out.versions.first()` idiom used throughout the subworkflows. It is cosmetic and pipeline-wide.

None of these affects the result. What decides pass or fail is the checks `tests/run_e2e.sh` prints at the end.

## What it does not cover

Protein annotation. That stage runs InterProScan, whose reference data is tens of gigabytes, so it cannot be part of a fast test. The profile skips it by supplying `representative_proteins_annotations`, which is the pipeline's own mechanism for providing that stage's output from elsewhere.

## The mock data

Three samples drawn from a shared three-genome synthetic community, each genome carrying twelve planted open reading frames. Committed under `assets/test_e2e` at roughly 0.3 MB in total. Regenerate with:

```bash
python tests/make_mock_metagenome.py --outdir assets/test_e2e
```

The generator is deterministic. It plants real ORFs rather than using random sequence, because Prodigal must find coding sequences or the clustering and abundance stages have no input, and it simulates reads from defined source sequences rather than at random, because MEGAHIT needs genuine overlap structure.

## Three repository defaults the genes route needs

Three settings are correctness-critical for the genes route and used to live only in `conf/metagear/*.config`, which `metagear-tools` merges at invocation but this repository's own `nextflow.config` does not load. A direct `nextflow run .` therefore failed, with errors that named none of the causes. Each now carries a repository default equal to the value the wrapper already sets, in `conf/modules.config` and `nextflow.config`:

- `VAMB_CONCATENATE_FASTA` needs `-m 10 --keepnames`. `bin/concatenate_fasta.py` otherwise applies a 2,000 bp floor that discards every called gene, and the run fails later at `Cannot seek at the end of file all.genes.fna`.
- `ABUNDANCE:COVERM_MAKE` needs `-p bwa-mem`. CoverM otherwise defaults to minimap2 and demands a FASTA beside the bwa index, which nf-core's `bwa/index` never writes there.
- `files_batch_size` needs a value. `abundance.nf` calls `bams.collate(params.files_batch_size)`, and a null argument raises `Unknown method invocation 'collate' on ArrayBag type` rather than behaving as unset.

Production behaviour is unchanged: the wrapper's configuration is merged last and continues to win. These defaults only decide what happens when nothing else sets them.

### Why the profile also sets `validation.ignoreParams`

`representative_proteins_annotations` and `amrfinder_db` are supplied by the wrapper and are deliberately declared in neither `nextflow.config` nor `nextflow_schema.json`, because nf-core's `schema_params` lint test fails on any parameter present in one of those files but not the other. The consequence is that nf-schema reports them as invalid input on a direct run, which looks alarming in the output of a test whose job is to reassure. The profile therefore lists those two names in `validation.ignoreParams`, scoped to the profile rather than pipeline-wide.

This suppresses only those two names. Samplesheet validation still runs: a samplesheet with an invalid `biome` is still rejected with the schema's own enum message, which was checked rather than assumed.
