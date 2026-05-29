# genes

De novo gene-centric analysis: assemble reads, call genes, build a non-redundant gene catalog, quantify it back across all samples, and group genes into metagenomic species pangenomes (MSPs). Use this when you need a sample-specific gene catalog or want to characterize organisms that are not in the MetaPhlAn reference.

## What it does

1. **Assembly** — MEGAHIT assembles each sample into contigs.
2. **Gene calling** — Prodigal (meta mode) predicts ORFs on the contigs; results are filtered to retain confidently-called genes.
3. **Gene clustering** — MMseqs2 (`--min-seq-id 0.95 -c 0.90`, alignment-mode 3) clusters all per-sample genes into a non-redundant catalog of representative genes.
4. **Protein translation** — representative genes are translated to proteins.
5. **Protein annotation** — AMRFinderPlus (`--plus`) annotates representative proteins against the AMR/virulence reference.
6. **Abundance quantification** — reads from each sample are aligned to the representative gene catalog with BWA, then CoverM produces three tables: raw counts, RPKM, and TPM (with `--min-read-percent-identity 95 --min-read-aligned-percent 75`).
7. **Optional MetaPhlAn profiling** — if `--metaphlan_profiles false` (the default), MetaPhlAn runs to produce per-sample taxonomic profiles used by MSP analysis. If a path is supplied instead, MetaPhlAn is skipped and those profiles are reused.
8. **MSP analysis** — MSPminer groups co-abundant genes into metagenomic species pangenomes; results are linked to GTDB-Tk taxonomy.

Each intermediate step can be skipped by supplying its output directly (see the "skip-mode" parameters below). This makes re-runs cheap when you want to vary only the last few steps.

## Inputs

- `--input` — samplesheet of QC'd reads.
- Databases: `--metaphlan_db`, `--gtdb_tk_db`, `--amrfinder_db` (from [download_databases](download_databases.md)).

## Parameters

| Parameter                                         | Type            | Default      | Controls                                                                                                     |
| ------------------------------------------------- | --------------- | ------------ | ------------------------------------------------------------------------------------------------------------ |
| `--input`                                         | path            | _(required)_ | Samplesheet of clean reads.                                                                                  |
| `--outdir`                                        | path            | _(required)_ | Result directory.                                                                                            |
| `--metaphlan_db`                                  | path            | —            | MetaPhlAn 4 database (used for MSP taxonomy).                                                                |
| `--gtdb_tk_db`                                    | path            | —            | GTDB-Tk reference for taxonomy.                                                                              |
| `--amrfinder_db`                                  | path            | —            | AMRFinderPlus database for protein annotation.                                                               |
| `--metaphlan_profiles`                            | path \| `false` | `false`      | If a path: directory of pre-computed MetaPhlAn profiles, MetaPhlAn skipped. If `false`: run MetaPhlAn fresh. |
| `--contigs_dir`                                   | path            | —            | **Skip assembly.** Directory containing `<sample>.contigs.fa.gz`.                                            |
| `--genes_dir`                                     | path            | —            | **Skip gene calling.** Directory containing `<sample>.all.genes.filtered.fasta`.                             |
| `--representative_genes`                          | path            | —            | **Skip clustering.** Use this representative gene catalog.                                                   |
| `--representative_proteins`                       | path            | —            | **Skip protein translation.** Use these representative proteins.                                             |
| `--representative_proteins_annotations`           | path            | —            | **Skip annotation.** Use these AMRFinderPlus annotations.                                                    |
| `--representative_genes_tpm` / `_rpkm` / `_count` | path            | —            | **Skip abundance** (all three must be set).                                                                  |

## Output

| Path (relative to `--outdir`)                                  | Content                                                                              |
| -------------------------------------------------------------- | ------------------------------------------------------------------------------------ |
| `assemblies/contigs/<sample>.contigs.fa.gz`                    | Per-sample assembled contigs (MEGAHIT).                                              |
| `assemblies/assembly_graphs/<sample>.k119.fastg.gz`            | Per-sample assembly de Bruijn graph (MEGAHIT).                                       |
| `prodigal/raw/`, `prodigal/`                                   | Raw and filtered Prodigal gene calls per sample (`*.all.genes.filtered.fasta`).      |
| `catalogs/raw/`                                                | Pre-clustering concatenated gene/protein sets (VAMB concatenation, translated proteins). |
| `catalogs/genes/all.genes.representative.fa.gz`                | **Representative gene catalog** (MMseqs2 cluster representatives).                   |
| `catalogs/genes/all.genes.clusters.tsv`                        | Gene cluster membership.                                                             |
| `catalogs/proteins/all.proteins.representative.fa.gz`          | **Representative protein catalog** (translation + MMseqs2 clustering).               |
| `catalogs/proteins/all.proteins.clusters.tsv`                  | Protein cluster membership.                                                          |
| `annotations/amrfinder/`                                       | AMRFinderPlus annotations on representative proteins.                                |
| `annotations/interproscan/`                                    | InterProScan annotations + FunctionalGroup parse.                                    |
| `abundance/all.genes/bwa_index/`                               | BWA-MEM index for the gene catalog.                                                  |
| `abundance/all.genes/bams/<sample>.bam`                        | Per-sample BAMs (reads mapped against the gene catalog).                             |
| `abundance/all.genes/per_batch/`                               | Per-batch CoverM contig tables (provenance for the merged matrices).                 |
| `abundance/all.genes/all.genes.count.tsv`                      | **Gene-by-sample raw read count matrix.**                                            |
| `abundance/all.genes/all.genes.rpkm.tsv`                       | **Gene-by-sample RPKM matrix.**                                                      |
| `abundance/all.genes/all.genes.tpm.tsv`                        | **Gene-by-sample TPM matrix.**                                                       |
| `metaphlan/individual_profiles/<sample>_microbial_profile.txt` | Per-sample MetaPhlAn profiles (only when MetaPhlAn runs).                            |
| `metaphlan/merged_microbial_profiles.txt`                      | Merged MetaPhlAn matrix (only when MetaPhlAn runs).                                  |
| `msp/`                                                         | MSPminer output: MSP definitions, pangenome sequences, abundance, MetaPhlAn linkage. |
| `pipeline_info/gene_analysis_multiqc_report.html`              | Consolidated MultiQC.                                                                |

The bolded rows are the typical analytical deliverables: the representative gene catalog, its annotations, the three abundance matrices, and the MSP set.

## Example

End-to-end run:

```bash
nextflow run schirmer-lab/metagear -profile docker \
  --workflow genes \
  --input clean.csv \
  --outdir genes/ \
  --metaphlan_db /data/metagear/metaphlan \
  --gtdb_tk_db /data/metagear/gtdb_tk \
  --amrfinder_db /data/metagear/amrfinder
```

Re-running only the MSP stage after upstream steps have completed:

```bash
nextflow run schirmer-lab/metagear -profile docker \
  --workflow genes \
  --input clean.csv \
  --outdir genes/ \
  --metaphlan_profiles previous_run/metaphlan/individual_profiles/ \
  --representative_genes previous_run/catalogs/genes/all.genes.representative.fa.gz \
  --representative_proteins previous_run/catalogs/proteins/all.proteins.representative.fa.gz \
  --representative_proteins_annotations previous_run/annotations/interproscan/all.proteins.FG_IPS_Pfam.tsv \
  --representative_genes_tpm previous_run/abundance/all.genes/all.genes.tpm.tsv \
  --representative_genes_rpkm previous_run/abundance/all.genes/all.genes.rpkm.tsv \
  --representative_genes_count previous_run/abundance/all.genes/all.genes.count.tsv \
  --metaphlan_db /data/metagear/metaphlan \
  --gtdb_tk_db /data/metagear/gtdb_tk \
  --amrfinder_db /data/metagear/amrfinder
```

## Notes

- **Compute footprint.** Assembly and clustering dominate; expect tens to hundreds of CPU-hours per sample, depending on diversity.
- **What "representative" means.** The MMseqs2 catalog clusters genes at ≥95% identity and ≥90% coverage. Representatives are the longest sequence per cluster — this is your effective gene unit for downstream stats.
- **Abundance filtering.** CoverM uses `--min-read-percent-identity 95 --min-read-aligned-percent 75` for counts and additionally `--min-covered-fraction 20` for RPKM/TPM. These thresholds are deliberately strict to limit cross-mapping between closely related genes.
- **Skip-mode design.** The optional inputs make re-runs incremental: if your catalog is stable but you want to redo the abundance step against more samples, only re-run from step 6 by passing the prior representatives.
- **MSPminer needs MetaPhlAn.** If you skip MetaPhlAn (`--metaphlan_profiles <dir>`), the MSP step uses the supplied profiles to link MSPs to taxonomy.
