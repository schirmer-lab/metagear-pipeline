---
title: Contributing
markdownPlugin: checklist
---

# `schirmer-lab/metagear-pipeline`: Contributing guidelines

Hi there!
Thanks for taking an interest in improving schirmer-lab/metagear-pipeline.

This page describes the recommended nf-core way to contribute to both schirmer-lab/metagear-pipeline and nf-core pipelines in general, including:

- [General contribution guidelines](#general-contribution-guidelines): common procedures or guides across all nf-core pipelines.
- [Pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines): procedures or guides specific to the development conventions of schirmer-lab/metagear-pipeline.

## General contribution guidelines

### Contribution quick start

To contribute code to any nf-core pipeline:

- [ ] Ensure you have Nextflow, nf-core tools, and nf-test installed. See the [nf-core/tools repository](https://github.com/nf-core/tools) for instructions.
- [ ] Check whether a GitHub [issue](https://github.com/schirmer-lab/metagear-pipeline/issues) about your idea already exists. If an issue does not exist, create one so that others are aware you are working on it.
- [ ] [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [schirmer-lab/metagear-pipeline repository](https://github.com/schirmer-lab/metagear-pipeline) to your GitHub account.
- [ ] Create a branch on your forked repository and make your changes following [pipeline conventions](#pipeline-contribution-conventions) (if applicable).
- [ ] To fix major bugs, name your branch `patch` and follow the [patch release](#patch-release) process.
- [ ] Update relevant documentation within the `docs/` folder, use nf-core/tools to update `nextflow_schema.json`, and update `CITATIONS.md`.
- [ ] Run and/or update tests. See [Testing](#testing) for more information.
- [ ] [Lint](#lint-tests) your code with nf-core/tools.
- [ ] Submit a pull request (PR) against the `dev` branch and request a review.

If you are not used to this workflow with Git, see the [GitHub documentation](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or [Git resources](https://try.github.io/) for more information.

### GitHub Codespaces

You can contribute to schirmer-lab/metagear-pipeline without installing a local development environment on your machine by using [GitHub Codespaces](https://github.com/codespaces).

[GitHub Codespaces](https://github.com/codespaces) is an online developer environment that runs in your browser, complete with VS Code and a terminal.
Most nf-core repositories include a devcontainer configuration, which creates a GitHub Codespaces environment specifically for Nextflow development.
The environment includes pre-installed nf-core tools, Nextflow, and a few other helpful utilities via a Docker container.

To get started, open the repository in [Codespaces](https://github.com/schirmer-lab/metagear-pipeline/codespaces).

### Testing

Once you have made your changes, run the pipeline with nf-test to test them locally.
For additional information, use the `--verbose` flag to view the Nextflow console log output.

```bash
nf-test test --tag test --profile +docker --verbose
```

If you have added new functionality, ensure you update the test assertions in the `.nf.test` files in the `tests/` directory.
Update the snapshots with the following command:

```bash
nf-test test --tag test --profile +docker --verbose --update-snapshots
```

When you create a pull request with changes, GitHub Actions will run automatic tests.
Pull requests are typically reviewed when these tests are passing.

Two types of tests are typically run:

#### Lint tests

nf-core has a [set of guidelines](https://nf-co.re/docs/specifications/overview) which all pipelines must follow.
To enforce these, run linting with nf-core/tools:

```bash
nf-core pipelines lint <pipeline_directory>
```

If you encounter failures or warnings, follow the linked documentation printed to screen.
For more information about linting tests, see [nf-core/tools API documentation](https://nf-co.re/docs/nf-core-tools/api_reference/latest/pipeline_lint_tests/actions_awsfulltest).

#### Pipeline tests

Each nf-core pipeline should be set up with a minimal set of test data.
GitHub Actions runs the pipeline on this data to ensure it runs through and exits successfully.
If there are any failures then the automated tests fail.
These tests are run with the latest available version of Nextflow and the minimum required version specified in the pipeline code.

### Patch release

> [!WARNING]
> Only in the unlikely event of a release that contains a critical bug.

- [ ] Create a new branch `patch` on your fork based on `upstream/main` or `upstream/master`.
- [ ] Fix the bug and use nf-core/tools to bump the version to the next semantic version, for example, `1.2.3` → `1.2.4`.
- [ ] Open a Pull Request from `patch` directly to `main`/`master` with the changes.

### Pipeline contribution conventions

nf-core semi-standardises how you write code and other contributions to make the schirmer-lab/metagear-pipeline code and processing logic more understandable for new contributors and to ensure quality.

#### Add a new pipeline step

To contribute a new step to the pipeline, follow the general nf-core coding procedure.
Please also refer to the [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines):

- [ ] Define the corresponding [input channel](#channel-naming-schemes) into your new process from the expected previous process channel.
- [ ] Install a module with nf-core/tools, or write a local module (see [default processes resource requirements](#default-processes-resource-requirements)), and add it to the target `<workflow>.nf`.
- [ ] Define the output channel if needed. Mix the version output channel into `ch_versions` and relevant files into `ch_multiqc`.
- [ ] Add new or updated parameters to `nextflow.config` with a [default value](#default-parameter-values).
- [ ] Add new or updated parameters and relevant help text to `nextflow_schema.json` with [nf-core/tools](#default-parameter-values).
- [ ] Add validation for relevant parameters to the pipeline utilisation section of `utils_nfcore_\_pipeline/main.nf` subworkflow.
- [ ] Perform local tests to validate that the new code works as expected.
  - [ ] If applicable, add a new test in the `tests` directory.
- [ ] Update `usage.md`, `output.md`, and `citation.md` as appropriate.
- [ ] [Lint](#lint-tests) the code with nf-core/tools.
- [ ] Update any diagrams or pipeline images as necessary.
- [ ] Update MultiQC config `assets/multiqc_config.yml` so relevant suffixes, file name cleanup, and module plots are in the appropriate order.
- [ ] If applicable, create a [MultiQC](https://seqera.io/multiqc/) module.
- [ ] Add a description of the output files and, if relevant, images from the MultiQC report to `docs/output.md`.

To update the minimum required Nextflow version, see the [Nextflow version bumping](#nextflow-version-bumping) section below. For more information about pipeline contributions, see [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines).

#### Channel naming schemes

Use the following naming schemes for channels to make the channel flow easier to understand:

- Initial process channel: `ch_output_from_<process>`
- Intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

#### Default parameter values

Parameters should be initialised and defined with default values within the `params` scope in `nextflow.config`.
They should also be documented in the pipeline JSON schema.

To update `nextflow_schema.json`, run:

```bash
nf-core pipelines schema build
```

The schema builder interface that loads in your browser should automatically update the defaults in the parameter documentation.

#### Default processes resource requirements

If you write a local module, specify a default set of resource requirements for the process.

Sensible defaults for process resource requirements (CPUs, memory, time) should be defined in `conf/base.config`.
Specify these with generic `withLabel:` selectors, so they can be shared across multiple processes and steps of the pipeline.

nf-core provides a set of standard labels that you should follow where possible, as seen in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/conf/base.config).
These labels define resource defaults for single-core processes, modules that require a GPU, and different levels of multi-core configurations with increasing memory requirements.

Values assigned within these labels can be dynamically passed to a tool using the the `${task.cpus}` and `${task.memory}` Nextflow variables in the `script:` block of a module (see an example in the [modules repository](https://github.com/nf-core/modules/blob/bd1b6a40f55933d94b8c9ca94ec8c1ea0eaf4b82/modules/nf-core/samtools/bam2fq/main.nf#L30)).

#### Nextflow version bumping

If you use a new feature from core Nextflow, bump the minimum required Nextflow version in the pipeline with:

```bash
nf-core pipelines bump-version --nextflow . <min_nf_version>
```

#### Images and figures guidelines

If you update images or graphics, follow the nf-core [style guidelines](https://nf-co.re/docs/community/brand/workflow-schematics).

## Pipeline specific contribution guidelines

A few conventions in this repository are not obvious from the code alone. Most of them exist because MetaGEAR is driven by an external wrapper, and several have silent failure modes — the pipeline runs, exits zero, and quietly does the wrong thing.

### Branching

The default branch is `master`; feature work targets `dev`. Feature pull requests into `dev` are squash-merged. Patch releases branch off `master` directly.

### Module include style

Includes use legacy `$projectDir` GString interpolation:

```groovy
include { FOO } from "$projectDir/path/to/module"
```

This is why CI sets `NXF_SYNTAX_PARSER=v1` and why the `nextflow-lint` pre-commit hook is disabled. **Do not migrate individual files to relative-path includes.** The move to v2 syntax is a coordinated change whose full scope is documented in `.pre-commit-config.yaml`; a piecemeal migration breaks the parser selection for everyone.

### Adding a parameter

Parameters reach the pipeline from two different places, and putting one in the wrong file breaks CI.

Parameters supplied by the wrapper — per-workflow inputs and database paths such as `contigs_dir`, `bins_dir`, `genomad_db` — are **not** declared in `nextflow.config` and **not** added to `nextflow_schema.json`. They arrive as `--param value` CLI flags and are documented in `workflow_definitions.json`. An undeclared parameter reads as `null`, which is exactly the default every consumer relies on (`if ( params.contigs_dir )`).

Only genuinely pipeline-level parameters belong in `nextflow.config`, and anything you add there **must** also be added to `nextflow_schema.json`. nf-core's `schema_params` lint test fails on any parameter present in one file but absent from the other, in either direction.

`workflow_definitions.json` is the source of truth for the user-facing CLI surface and is consumed by the wrapper. Keep it in sync when adding or renaming workflow parameters.

### Resource requests

Edit `conf/resources.yaml`, then regenerate `conf/resources.config` with the wrapper's generator:

```bash
bash <metagear-tools>/lib/build_resources_config.sh conf/resources.yaml conf/resources.config
```

Never hand-edit `conf/resources.config` — it carries an AUTO-GENERATED banner and your change will be lost on the next regeneration.

> [!WARNING]
> Each entry in `conf/resources.yaml` **must** open and close its `{ ... }` flow mapping on a single physical line. The generator parses the file line by line; an entry split across lines is silently skipped, with no warning and a zero exit status, and the process falls back to bare label defaults. This has already caused one real regression. `conf/resources.yaml` is listed in `.prettierignore` for exactly this reason — do not remove it, and do not run a formatter over that file.

### Per-workflow configs

The files under `conf/metagear/` are **not** loaded by this repository's `nextflow.config`. The wrapper merges them at invocation time. Keep each one self-contained, and remember that when you run `nextflow run .` directly they are inactive — a `withName` override you add there will not apply to a direct invocation.

Because the wrapper merges _all_ of them regardless of the selected workflow, a `ext.args` GString referencing a parameter is evaluated for every workflow. Such parameters belong in `conf/metagear/common.config` rather than `nextflow.config`, so the declaration and the reference share a load cycle. `drep_batch_ani` and `iphop_batch_size` follow this pattern.

### Publishing outputs

The results tree is additive across workflows — see [output.md](output.md) for the layout and naming conventions. Prefer `publishDir` and `saveAs` changes over renaming a module's emitted filenames: publish-time renames do not alter the process hash, so existing runs still resume cleanly.

### Before opening a pull request

Run the same three checks CI runs:

```bash
nf-test test --profile +docker            # add --update-snapshots after intentional output changes
nf-core pipelines lint
prek run --all-files                      # or: pre-commit run --all-files
```

If you changed `manifest.name`, expect knock-on lint failures: nf-core derives logo filenames, MultiQC `report_section_order` section IDs, and issue-template links from it, and the nf-test version snapshot embeds it.

New dependencies in `bin/` scripts must already exist in the container that runs them — these scripts execute inside the per-process container, not on the host.
