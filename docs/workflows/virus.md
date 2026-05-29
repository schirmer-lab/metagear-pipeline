# virus

End-to-end viral and plasmid analysis: detect viral/plasmid contigs from a metagenome assembly, cluster them into non-redundant catalogs, annotate them, predict their bacterial hosts, and identify auxiliary metabolic genes (AMGs). Use this when your question is about phages, plasmids, mobile genetic elements, or the metabolic functions viruses confer on their hosts.

## What it does

1. **Assembly** — MEGAHIT, same as in [genes](genes.md). Skippable with `--contigs_dir`.
2. **Two-pass viral detection.** geNomad and CheckV run in two passes:
   - **Pass 1** — `GENOMAD_PASS1` (`--cleanup --conservative --enable-score-calibration`) classifies viral and plasmid contigs; `CHECKV_PASS1` (`--remove_tmp`) trims host flanks and reports completeness.
   - **Pass 2** — runs again on the trimmed proviruses to recover viral signal that was masked by host sequence. Outputs are merged and filtered through ICTV taxonomy.
3. **Viral and plasmid clustering** — MMseqs2 clusters viral contigs (`--min-seq-id 0.95 -c 0.85`, easy-linclust) and plasmid contigs (`--min-seq-id 0.95 -c 0.90`) into non-redundant catalogs.
4. **Gene calling and clustering** — Prodigal calls genes on all contigs; MMseqs2 clusters them; representative genes are translated to proteins, which are themselves clustered (`--min-seq-id 0.90 -c 0.80`).
5. **Protein annotation** — AMRFinderPlus (`--plus`) on the representative proteins.
6. **Viral gene representatives** — for each viral and plasmid contig, the gene IDs are pulled from the full clustered set and emitted as a separate "viral genes" and "plasmid genes" catalog, with annotations preserved.
7. **Abundance quantification** — CoverM aligns reads against the viral, plasmid, and gene catalogs to produce count/RPKM/TPM tables (`--min-read-percent-identity 95 --min-read-aligned-percent 75`, `--min-covered-fraction 20`).
8. **Viral annotation** — Pharokka (`--mmseqs2_only`) annotates phage proteins; VirSorter2 (`--prep-for-dramv`, dsDNAphage/ssDNA/NCLDV/lavidaviridae, length ≥ 1500) preprocesses for DRAM-v; DRAM-v (`--skip_trnascan --min_contig_size 1500`) identifies AMGs; iPHoP predicts bacterial hosts at genus and genome level.
9. **AMG post-processing** — AMGs are joined with the abundance tables to produce AMG-specific count/RPKM/TPM matrices.

The two-pass geNomad + CheckV is the central novelty: it recovers viral sequences embedded in host genomes (proviruses) that a single pass would miss.

## Inputs

- `--input` — samplesheet of QC'd reads.
- Databases: `--genomad_db`, `--checkv_db`, `--virsorter2_db`, `--pharokka_db`, `--dram_db`, `--iphop_db`, `--amrfinder_db` (from [download_databases](download_databases.md) with `--databases virus` or `all`).

## Parameters

| Parameter         | Type | Default      | Controls                                                                 |
| ----------------- | ---- | ------------ | ------------------------------------------------------------------------ |
| `--input`         | path | _(required)_ | Samplesheet of clean reads.                                              |
| `--outdir`        | path | _(required)_ | Result directory.                                                        |
| `--genomad_db`    | path | —            | geNomad database.                                                        |
| `--checkv_db`     | path | —            | CheckV database.                                                         |
| `--virsorter2_db` | path | —            | VirSorter2 database.                                                     |
| `--pharokka_db`   | path | —            | Pharokka database.                                                       |
| `--dram_db`       | path | —            | DRAM-v database.                                                         |
| `--iphop_db`      | path | —            | iPHoP database.                                                          |
| `--amrfinder_db`  | path | —            | AMRFinderPlus database.                                                  |
| `--contigs_dir`   | path | —            | **Skip assembly.** Directory of `<sample>.contigs.fa.gz`.                |
| `--genes_dir`     | path | —            | **Skip gene calling.** Directory of `<sample>.all.genes.filtered.fasta`. |

## Output

| Path (relative to `--outdir`)                         | Content                                                                |
| ----------------------------------------------------- | ---------------------------------------------------------------------- |
| `assemblies/contigs/<sample>.contigs.fa.gz`           | Per-sample contigs (MEGAHIT).                                          |
| `assemblies/assembly_graphs/<sample>.k119.fastg.gz`   | Per-sample assembly de Bruijn graph (MEGAHIT).                         |
| `virus/genomad/virus/<sample>/`                       | geNomad pass-1 output (viral classification).                          |
| `virus/checkv/virus/<sample>/virus.*`                 | CheckV pass-1 output (completeness/contamination).                     |
| `virus/genomad/provirus/<sample>/`                    | geNomad pass-2 (on trimmed proviruses).                                |
| `virus/checkv/provirus/<sample>/provirus.*`           | CheckV pass-2.                                                         |
| `virus/per_sample/<sample>/`                          | Per-sample viral detection intermediates: `virus.summary[.filtered].tsv`, `virus.ids.txt`, `virus.filtered.fna.gz`, plus the plasmid counterparts. |
| `catalogs/contigs/virus.contigs.representative.fa.gz` | **Non-redundant viral catalog (vOTUs).**                               |
| `catalogs/contigs/virus.contigs.clusters.tsv`         | Viral contig cluster membership.                                       |
| `catalogs/contigs/plasmid.contigs.representative.fa.gz` | **Non-redundant plasmid catalog.**                                   |
| `catalogs/contigs/plasmid.contigs.clusters.tsv`       | Plasmid contig cluster membership.                                     |
| `catalogs/genes/all.genes.representative.fa.gz`       | Cohort gene catalog (same form as `genes`).                    |
| `catalogs/genes/all.genes.clusters.tsv`               | Cohort gene cluster membership.                                        |
| `catalogs/genes/virus.genes.representative.fa.gz`     | Virus-only gene catalog.                                               |
| `catalogs/genes/virus.genes.clusters.tsv`             | Virus-only gene cluster membership.                                    |
| `catalogs/genes/plasmid.genes.representative.fa.gz`   | Plasmid-only gene catalog.                                             |
| `catalogs/genes/plasmid.genes.clusters.tsv`           | Plasmid-only gene cluster membership.                                  |
| `catalogs/proteins/all.proteins.representative.fa.gz` | Cohort protein catalog.                                                |
| `catalogs/proteins/all.proteins.clusters.tsv`         | Cohort protein cluster membership.                                     |
| `catalogs/proteins/virus.proteins.representative.fa.gz`   | Virus-only protein catalog.                                        |
| `catalogs/proteins/plasmid.proteins.representative.fa.gz` | Plasmid-only protein catalog.                                      |
| `catalogs/raw/`                                       | Pre-clustering concatenated FASTAs (intermediates).                    |
| `clusters/all.genes.clusters.annotated.tsv`           | Cohort gene-cluster annotations with per-member class labels.          |
| `clusters/all.genes.clusters.aggregated.tsv`          | Per-representative rollup of cluster annotations.                      |
| `virus/clusters/<scope>.genes.clusters.annotated.tsv` | Virus/plasmid-scoped gene cluster annotations.                         |
| `virus/clusters/<scope>.genes.representative_ids.txt` | Reps for clusters touching `<scope>` (virus or plasmid).               |
| `virus/clusters/<scope>-exclusive.genes.representative_ids.txt` | Reps for clusters whose members are entirely `<scope>`.            |
| `virus/pharokka/`                                     | Pharokka annotations (`*.gff`, `*.fna`, `*.faa`).                      |
| `virus/virsorter2/`                                   | VirSorter2 output prepared for DRAM-v.                                 |
| `virus/dramv/parts/`, `virus/dramv/`                  | DRAM-v annotations and AMG calls.                                      |
| `virus/iphop/parts/`                                  | iPHoP host predictions (genus and genome).                             |
| `abundance/virus/virus.{count,rpkm,tpm}.tsv`          | **Viral abundance matrices** (vOTU × sample).                          |
| `abundance/virus/virus.amg.{count,rpkm,tpm}.tsv`      | **AMG abundance matrices** (per-sample × AMG counts).                  |
| `abundance/plasmid/plasmid.{count,rpkm,tpm}.tsv`      | **Plasmid abundance matrices.**                                        |
| `abundance/all.genes/all.genes.{count,rpkm,tpm}.tsv`  | Gene-level abundance matrices.                                         |
| `abundance/<class>/bams/SAMPLE-N.bam`                 | Per-sample BAMs (reads mapped against the catalog).                    |
| `abundance/<class>/bwa_index/`                        | BWA-MEM index for the catalog.                                         |
| `abundance/<class>/per_batch/batch_NNN.<metric>.tsv`  | Per-batch coverm-contig matrices — provenance for the merged tables.   |
| `virus/detection/virus.calls.tsv`                     | **All per-sample viral calls** (geNomad + CheckV, cohort-concatenated). |
| `virus/detection/virus.calls.filtered.tsv`            | **FDR-passing per-sample viral calls.**                                |
| `virus/detection/plasmid.calls.tsv`                   | All per-sample plasmid calls.                                          |
| `virus/detection/plasmid.calls.filtered.tsv`          | **FDR-passing per-sample plasmid calls.**                              |
| `virus/lifestyle/virus.amgs.tsv`                      | **AMG identifications on the vOTU catalog** (DRAM-V).                  |
| `virus/lifestyle/virus.amg_to_catalog.tsv`            | **AMG → gene catalog cross-reference** (primary deliverable).          |
| `virus/lifestyle/virus.host.genome.tsv`               | iPHoP host predictions, genome-level (per vOTU).                       |
| `virus/lifestyle/virus.host.genus.tsv`                | iPHoP host predictions, genus-level (per vOTU).                        |
| `virus/dramv/`                                        | Supporting DRAM-V artefacts (`amg_filtered.{fna,faa}`, search results, ID lists, per-batch parts). |
| `pipeline_info/viral_analysis_multiqc_report.html`    | Consolidated MultiQC.                                                  |

The bolded rows are the primary deliverables for most downstream analyses.

## Example

```bash
nextflow run schirmer-lab/metagear -profile docker \
  --workflow virus \
  --input clean.csv \
  --outdir viruses/ \
  --genomad_db /data/metagear/genomad \
  --checkv_db /data/metagear/checkv \
  --virsorter2_db /data/metagear/virsorter2 \
  --pharokka_db /data/metagear/pharokka \
  --dram_db /data/metagear/dram \
  --iphop_db /data/metagear/iphop \
  --amrfinder_db /data/metagear/amrfinder
```

Sharing an assembly with a prior `genes` run:

```bash
nextflow run schirmer-lab/metagear -profile docker \
  --workflow virus \
  --input clean.csv \
  --outdir viruses/ \
  --contigs_dir prior_run/assembly/ \
  --genes_dir prior_run/prodigal/ \
  ...
```

## Notes

- **Why two passes.** A single geNomad + CheckV pass catches free viral contigs but misses proviruses still embedded in host sequence. Pass 2 runs on CheckV-trimmed proviruses to recover that signal. Expect noticeably more recovered viral sequence than a single-pass pipeline.
- **VirSorter2 is used only for DRAM-v prep**, not for primary detection. Detection is geNomad + CheckV. The VirSorter2 invocation (`--include-groups dsDNAphage,ssDNA,NCLDV,lavidaviridae --min-length 1500 --min-score 0.5`) is tuned for that prep role.
- **Compute and memory.** This is the heaviest workflow in the suite — DRAM-v and iPHoP are both demanding. Expect substantial wall time even on small studies; consider running it on a subset of samples first.
- **AMGs are auxiliary metabolic genes** — phage-encoded genes that augment host metabolism. The AMG post-processing step joins DRAM-v's AMG calls with the gene abundance matrices, producing AMG-specific TPM/RPKM/count tables suitable for cross-sample analysis.
- **Plasmids ride along.** geNomad classifies plasmids alongside viruses; the plasmid catalog and its abundance tables are emitted in parallel with the viral set. If you only care about viruses, ignore the plasmid outputs.
