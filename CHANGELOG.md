# schirmer-lab/metagear-pipeline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versions are date-based (`YY.MM`, or `YY.MM.N` for a fix within the same month)
and are released in lockstep with [metagear-tools](https://github.com/schirmer-lab/metagear-tools).

## v26.09.1 - [2026-09-01]

Documentation only. No pipeline code changes; results are identical to 26.09.

### `Added`

- Reference pages for `classification`, `mag`, `msp` and `structures`, the four workflows that
  shipped in 26.09 without one. `docs/workflows/index.md` now links all ten.

### `Fixed`

- The README listed four workflows and named `gene_analysis`, which 26.09 renamed to `genes`.
- The README quick-start omitted `--workflow`. Without it a run validates the samplesheet and
  exits without analysing anything, so the example appeared to work while doing nothing.

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
