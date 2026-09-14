# schirmer-lab/metagear-pipeline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versions are date-based (`YY.MM`, or `YY.MM.N` for a fix within the same month)
and are released in lockstep with [metagear-tools](https://github.com/schirmer-lab/metagear-tools).

## v26.10dev - [unreleased]

### `Added`

- `easy_map`: map reads against a catalog the caller supplies. It assembles nothing and
  calls no genes, so the catalog can come from anywhere, including another pipeline. Takes
  `--input` and `--catalog`, plus `--catalog_label` to name the outputs, and publishes
  `abundance/<label>/<label>.{count,covered_bases,rpkm,tpm}.tsv` in the same shape as the
  other abundance classes.
- The abundance metrics are now configuration. `COVERM_CONTIG` and `COVERM_GENOME` take
  `ext.methods` (computed with `ext.args`) and `ext.methods2` (computed with `ext.args2`),
  one coverm call and one published table per metric, so adding or dropping a metric no
  longer touches the modules or `ABUNDANCE`.
- Breadth of coverage for contig-mode callers (`genes`, `virus`):
  `abundance/<class>/<class>.covered_bases.tsv`, the number of reference bases covered per
  feature per sample. It sits in `ext.methods` beside `count`, without
  `--min-covered-fraction`, so a zero means no covered base rather than a value below a
  threshold, which is what a detection criterion needs.

### `Changed`

- `trimmed_mean` is no longer computed. Both coverm modules produced it on every run and
  `ABUNDANCE` never consumed it, so it was never published. It can be restored for any
  caller by adding it to `ext.methods2`, and it will then be published like any other
  metric.
- Comments across `conf/` now carry only what the code cannot state. Removed were a
  published-layout tree that contradicted the block below it, descriptions of output
  files that `docs/output.md` owns, restatements of the rename each `saveAs` performs,
  superseded `ext.args` kept as history, and a claim that the raw gene classification
  table holds virus labels only, which it does not. Kept were tool quirks, thresholds
  with a reason, and resource numbers a failure justified.
- Removed an empty `process { }` block in `conf/metagear/genes.config` left behind when
  the MSP blocks moved out. It made Nextflow warn `Unknown config attribute 'process'`
  on every run.
- `COVERM_CONTIG_BATCH` is removed. Batching moved into `ABUNDANCE` itself, which splits
  each label's BAMs with `collate(params.files_batch_size)`, and the process had been an
  unused include since `d72ff74` with its own hardcoded batch size of 50.

## v26.09 - [2026-09-01]

First release of the integrated microbiome pipeline. This is a major version and
it is not backwards compatible with 1.0.1 — workflow names, the results tree and
several parameters have changed. The 1.x line remains available at
`schirmer-lab/metagear-pipeline-legacy`.

### `Added`

- `virus` — viral and plasmid detection, vOTU clustering (vclust, MIUViG criterion),
  and annotation via Pharokka, VirSorter2, DRAM-v, iPHoP and PhaTYP lifestyle calls.
- `classification` — viral/plasmid partitioning with geNomad, bacterial binning with
  SemiBin2 + MetaBAT2 refined by Binette and scored by CheckM2, and a per-contig
  classification table. Reads an optional `biome` samplesheet column to pick the
  SemiBin2 model.
- `mag` — cohort MAG catalog: dRep dereplication, GTDB-Tk taxonomy on the
  representatives, and MAG×sample abundance.
- `msp` — MetaSpecies Pangenomes via MSPminer co-abundance clustering, with GTDB-Tk
  taxonomy and a MetaPhlAn cross-walk.
- `structures` — protein structural-homology annotation through PHOLD
  (ProstT5 → Foldseek), including an offsite-GPU packaging mode.
- Presets that chain the workflows in one command: `profiles`, `genomes`, `microbiome`.
- `--reuse-outputs` support across workflows, so a later workflow discovers an
  earlier one's artifacts instead of recomputing them.

### `Changed`

- **Breaking** — `gene_analysis` is now `genes`, and its catalog is built with MMseqs2
  rather than CD-HIT, so gene and protein representatives differ from 1.x.
- **Breaking** — the results tree is reorganised: catalogs, abundance, annotations,
  assemblies and per-workflow directories replace the previous layout.
- `pipeline_info` outputs are prefixed with the workflow name, so several workflows
  can share one `--outdir` without their reports being indistinguishable.
- nf-core template updated to 4.0.2; Nextflow >= 25.10.4 is now required.

### `Fixed`

- `virus` re-derived the gene and protein catalogs even when they were supplied.
  Because MMseqs2 picks different cluster representatives each run, this silently
  invalidated the gene-cluster classification table and the MSP catalog.
- DAG construction aborted in `classification` and `mag` by reading `.out.versions`
  on processes that publish versions to a topic channel.
- DAG construction aborted in `genes` when `--representative_proteins` was supplied
  without `--representative_proteins_annotations`.
