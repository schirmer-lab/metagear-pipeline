# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Pipeline overview

`schirmer-lab/metagear-pipeline` is an nf-core-templated Nextflow DSL2 pipeline for shotgun metagenomics. It is a **multi-workflow pipeline** — `main.nf` always enters `workflows/metagear.nf`, which dispatches to one of several entry-point workflows based on `params.workflow`:

| `params.workflow` value                         | Entry subworkflow                                                                                                                   |
| ----------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `download_databases`                            | `workflows/setup.nf` (`SETUP`) — fetches Kneaddata, MetaPhlAn, HUMAnN, etc.                                                         |
| `qc_dna`, `qc_rna` (or anything starting `qc_`) | `subworkflows/local/common/quality_control.nf`                                                                                      |
| `microbial_profiles`                            | `subworkflows/local/microbiome/microbial_profiles.nf` (MetaPhlAn + HUMAnN)                                                          |
| `genes`                                         | `subworkflows/local/microbiome/genes.nf` (gene catalog: assembly → gene call → clustering → abundance → protein annotation)         |
| `virus`                                         | `subworkflows/local/microbiome/virus.nf`                                                                                            |
| `classification`                                | `subworkflows/local/microbiome/classification.nf` (assembly + viral/plasmid partition + bacterial binning + per-contig TSV)         |
| `mag`                                           | `subworkflows/local/microbiome/mag.nf` (cohort MAG catalog: dRep + GTDB-Tk + coverm-genome; builds on a prior `classification` run) |
| `msp`                                           | `subworkflows/local/pangenome/msp.nf` (MetaSpecies Pangenomes: MSPminer + GTDB-Tk + MetaPhlAn; builds on a prior `genes` run)       |

Every entry subworkflow follows a `*_INIT` + `*` pair: the `_INIT` validates inputs and constructs DB channels; the main workflow takes those channels and runs the analysis. `SUMMARY` then collects `versions` and `summary_data` channels for MultiQC + version reporting.

`workflow_definitions.json` is the source of truth for the user-facing CLI surface (workflows and their parameters); it is consumed by the external `metagear-tools` wrapper. Keep it in sync when adding or renaming workflow params.

## Common commands

```bash
# Run the test profile (uses conf/test.config; pins params.workflow = "test")
nextflow run . -profile test,docker --outdir results

# Run a specific workflow
nextflow run . -profile docker --workflow genes --input samplesheet.csv --outdir results

# nf-test — execute the pipeline test suite (matches the CI matrix)
nf-test test --tag test --profile +docker --verbose

# Update snapshots after intentional output changes
nf-test test --tag test --profile +docker --update-snapshots

# nf-core linting (CI runs this; --release on PRs into master/main)
nf-core pipelines lint

# Pre-commit hooks (prettier + whitespace fixers). Nextflow lint is intentionally disabled — see .pre-commit-config.yaml
pre-commit run --all-files
```

`nf-test` runs against `tests/default.nf.test`. The snapshot it asserts on is `pipeline_info/test_mqc_versions.yml` — the `test_` prefix comes from `params.workflow = "test"` pinned in `conf/test.config` (do not set this in the `.nf.test` file; `workflow` is a reserved identifier there).

## Configuration architecture

`nextflow.config` loads, in order: `conf/base.config` (default per-process resources + retry policy), `conf/resources.config` (resource tiers and per-process overrides), then an optional user override at `$INSTALL_DIR/resources.config` (or `~/.metagear/resources.config`), and finally `conf/modules.config` (publishDir defaults).

Important subsystems:

- **Resource model.** Edit `conf/resources.yaml` and regenerate `conf/resources.config` via the `metagear-tools` generator (`build_resources_config.sh`). Do not hand-edit `resources.config` — it has an "AUTO-GENERATED" banner. Memory/time get `* task.attempt` for retry escalation; `process.resourceLimits` (set from `params.max_cpus/max_memory/max_time` in `nextflow.config`) caps every process.
- **`max_cpus` / `max_memory` / `max_time` are intentional params** even though nf-core 3.x deprecated them. They are wired into `process.resourceLimits` and `executor.cpus/memory` so the external wrapper can override them per-user. See `.nf-core.yml` (the lint ignore-list documents the rationale).
- **`conf/metagear/*.config`** (per-workflow ext.args / publishDir overrides) are **not loaded by this repo's `nextflow.config`**. They are merged in by the external `metagear-tools` wrapper at invocation time. Keep them self-contained — when running `nextflow run .` directly they will not be active.

## Code layout conventions

- `modules/local/<tool>/main.nf` — pipeline-specific process definitions. `modules/nf-core/` is vendored from nf-core/modules and managed by `nf-core modules` (tracked in `modules.json`); do not hand-edit those.
- `subworkflows/local/common/` — building blocks reused across entry workflows (input check, assembly, gene call, clustering, abundance, protein annotation, summary).
- `subworkflows/local/{microbiome,virus,pangenome,setup}/` — domain-specific pipelines composed from `common/`.
- `bin/` — Python/shell helpers invoked from processes. They run inside the per-process container, so any new dependency must already be available in that container.
- Module include style is `include { FOO } from "$projectDir/path/to/module"` (legacy `$projectDir` GString interpolation). This is why CI sets `NXF_SYNTAX_PARSER=v1` and `.pre-commit-config.yaml` has the nextflow-lint hook disabled. Do not switch new code to relative-path includes piecemeal — the migration to v2 syntax is a coordinated change documented in the pre-commit config.

## Testing and CI

- `.github/workflows/nf-test.yml` shards tests across 7 shards and runs each on `docker` (always) plus `singularity` (only for PRs into `master`/`main`). The `conda` profile is currently excluded; re-enable only after verifying all module-level conda envs solve. Tests run on both Nextflow `25.10.4` and `latest-everything` (the latter is non-blocking).
- `.github/workflows/linting.yml` runs `nf-core pipelines lint` (with `--release` on PRs into master/main) and `prek` (pre-commit replacement) for prettier + whitespace.
- Minimum Nextflow version is pinned in `nextflow.config` (`!>=25.10.4`). Bump via `nf-core pipelines bump-version --nextflow . <ver>`.

## Branching

Default branch is `master`; feature work targets `dev`. Patch releases branch off `master` directly per `docs/CONTRIBUTING.md`.
