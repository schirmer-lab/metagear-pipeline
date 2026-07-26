# Workflows

The `schirmer-lab/metagear-pipeline` pipeline groups its work into ten entry-point workflows, selected at run time with the `--workflow` parameter:

| Workflow                                    | Purpose                                                                             | Input                   | Output                                                            | Cost           |
| ------------------------------------------- | ----------------------------------------------------------------------------------- | ----------------------- | ----------------------------------------------------------------- | -------------- |
| [download_databases](download_databases.md) | One-time install of all reference databases                                         | —                       | KneadData, MetaPhlAn4, HUMAnN3, GTDB-Tk, plus viral databases     | Disk + network |
| [qc_dna](qc_dna.md)                         | Adapter/quality trimming and host decontamination of DNA reads                      | Raw DNA FASTQ           | Clean paired reads + QC report                                    | Medium         |
| [qc_rna](qc_rna.md)                         | Same flow as `qc_dna`; intended for metatranscriptomic input                        | Raw RNA FASTQ           | Clean paired reads + QC report                                    | Medium         |
| [microbial_profiles](microbial_profiles.md) | Reference-based taxonomic and functional profiling                                  | Clean reads             | MetaPhlAn4 species table + HUMAnN3 gene-family and pathway tables | Medium–High    |
| [genes](genes.md)                           | De novo assembly, gene calling, gene catalog, MSP analysis                          | Clean reads             | Gene/protein representative catalogs, abundance matrices, MSPs    | High           |
| [virus](virus.md)                           | Viral and plasmid detection, clustering, annotation, host prediction, AMG discovery | Clean reads             | Viral and plasmid catalogs, AMGs, iPHoP host predictions          | Very high      |
| `classification`                            | Viral/plasmid partition, bacterial binning, per-contig classification               | Clean reads             | Per-contig classification TSV, per-sample MAG bins                | Very high      |
| `mag`                                       | Cohort MAG catalog — dRep, GTDB-Tk taxonomy, MAG×sample abundance                   | `classification` output | MAG catalog, GTDB-Tk lineages, abundance matrices                 | High           |
| `msp`                                       | MetaSpecies Pangenomes — MSPminer co-abundance clustering, GTDB-Tk, MetaPhlAn       | `genes` output          | MSP membership, taxonomy, MSP×sample abundance                    | High           |
| `structures`                                | Protein structural-homology annotation via PHOLD (ProstT5 → Foldseek)               | `genes`/`virus` output  | Per-representative structural annotations                         | High (GPU)     |

> [!NOTE]
> Dedicated pages for `classification`, `mag`, `msp`, and `structures` are still being written. Until they land, `metagear <workflow> --help` and `workflow_definitions.json` are the authoritative parameter references, and [output.md](../output.md) documents what each produces.

## Recommended order

Run `download_databases` once per machine, then quality-control the raw reads, then choose whichever downstream workflows match your scientific question. The four read-level analysis workflows (`microbial_profiles`, `genes`, `virus`, `classification`) are complementary and can run on the same QC'd reads.

The remaining three are second-stage workflows that read a previous run's outputs from disk rather than starting from reads: `mag` follows `classification`, `msp` follows `genes`, and `structures` follows either `genes` or `virus`. Point them at the same `--outdir` as the run they build on, or let the wrapper's `--reuse-outputs` discover the inputs automatically.

```bash
# 1. One-time database setup
nextflow run schirmer-lab/metagear-pipeline -profile docker \
  --workflow download_databases --outdir databases/

# 2. QC the raw reads
nextflow run schirmer-lab/metagear-pipeline -profile docker \
  --workflow qc_dna --input raw.csv --outdir qc/

# 3. Pick one or more analyses
nextflow run schirmer-lab/metagear-pipeline -profile docker \
  --workflow microbial_profiles --input clean.csv --outdir profiles/

nextflow run schirmer-lab/metagear-pipeline -profile docker \
  --workflow genes --input clean.csv --outdir genes/

nextflow run schirmer-lab/metagear-pipeline -profile docker \
  --workflow virus --input clean.csv --outdir viruses/
```

## Picking a workflow

- **Which organisms are present, and what functions can they perform?** → `microbial_profiles` (reference-based, no assembly, fast).
- **What genes are in this sample, including novel ones?** → `genes` (assembly-based, captures genes absent from reference databases, builds metagenomic species pangenomes).
- **What viruses, phages, plasmids, or auxiliary metabolic genes are present?** → `virus` (assembly-based, two-pass viral detection with geNomad + CheckV).
- **What is each contig, and which genomes can I recover?** → `classification`, then `mag` for a dereplicated cohort MAG catalog with GTDB-Tk taxonomy.
- **Which species-level gene groups co-vary across my cohort?** → `genes`, then `msp`.
- **What do the hypothetical proteins actually do?** → `structures`, which recovers function for proteins that sequence-homology annotation leaves unassigned.

## Input samplesheet

All workflows except `download_databases` take a CSV samplesheet via `--input`. Columns: `sample`, `fastq_1`, `fastq_2` (R2 optional for single-end), plus an optional `biome` column between `sample` and `fastq_1` used by `classification`. See [usage.md](../usage.md) for the full samplesheet contract.

## Wrapper CLI

The [`metagear-tools`](https://github.com/schirmer-lab/metagear-tools) CLI wraps these workflows so that day-to-day invocations look like `metagear qc_dna --input samples.csv`. The pipeline runs identically with or without the wrapper; the Nextflow invocations shown here are the canonical form.
