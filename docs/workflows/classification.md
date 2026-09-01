# classification

Decide what every assembled contig is, and recover bacterial genomes from the chromosomal fraction. Contigs are partitioned into viral, plasmid and chromosomal sets with geNomad, the chromosomal set is binned by two independent binners and refined into MAGs, and everything is reconciled into one per-contig table. Use this when you want genome-resolved structure for a cohort, or a per-contig answer to "is this bacterial, viral, plasmid, or unknown?".

## What it does

1. **Assembly** — MEGAHIT assembles each sample into contigs, unless `--contigs_dir` supplies them.
2. **Viral/plasmid partition** — geNomad runs in two passes with CheckV, splitting each assembly into viral, plasmid and chromosomal sequences. `--chromosome_dir` skips this entirely when a previous `virus` run already produced the chromosome partition.
3. **Read mapping** — each sample's reads are mapped back to its own chromosomal contigs with BWA, and CoverM produces the BAM used for coverage.
4. **Binning, twice** — MetaBAT2 (`--minContig 1500`) bins on contig depth, and SemiBin2 bins with a pretrained model chosen from the samplesheet's `biome` column.
5. **Bin refinement** — Binette reconciles the two bin sets and scores them with CheckM2, keeping bins at `--min_completeness 50` (MIMAG medium-quality and above).
6. **Eukaryote screen** — Tiara labels contigs as prokaryotic, eukaryotic, organellar or unknown, which is what flags host and dietary sequence that geNomad has no opinion about.
7. **Unbinned taxonomy (optional)** — with `--enable_contig_taxonomy`, MMseqs2 `easy-taxonomy` assigns lineage to contigs that ended up in no bin. Off by default.
8. **Reconciliation** — every evidence source is merged into one per-contig table, and the gene catalog is cross-walked onto it when `--gene_clusters_tsv` is supplied.

## Inputs

- `--input` — samplesheet of QC'd reads. May carry an optional `biome` column between `sample` and `fastq_1`.
- Databases: `--genomad_db`, `--checkv_db`, `--checkm2_db`, and `--mmseqs_taxonomy_db` only if contig taxonomy is enabled (from [download_databases](download_databases.md)).

## Parameters

| Parameter                  | Type | Default      | Controls                                                                                          |
| -------------------------- | ---- | ------------ | ------------------------------------------------------------------------------------------------- |
| `--input`                  | path | _(required)_ | Samplesheet of clean reads.                                                                       |
| `--outdir`                 | path | _(required)_ | Result directory.                                                                                 |
| `--genomad_db`             | path | —            | geNomad database.                                                                                 |
| `--checkv_db`              | path | —            | CheckV database.                                                                                  |
| `--checkm2_db`             | path | —            | CheckM2 database, used by Binette for bin scoring.                                                |
| `--contigs_dir`            | path | —            | **Skip assembly.** Directory of `<sample>.contigs.fa.gz`.                                         |
| `--chromosome_dir`         | path | —            | **Skip viral detection.** Directory of `<sample>.chromosome.fna.gz` from a prior run.             |
| `--viral_ids_dir`          | path | —            | Viral contig IDs from a prior `virus` run, so the per-contig table can still label virus contigs. |
| `--plasmid_ids_dir`        | path | —            | Plasmid contig IDs, same purpose.                                                                 |
| `--gene_clusters_tsv`      | path | —            | Gene clusters from a prior `genes` run; enables the gene-catalog cross-walk.                      |
| `--enable_contig_taxonomy` | bool | `false`      | Run MMseqs2 taxonomy on unbinned contigs. Expensive; see Notes.                                   |
| `--mmseqs_taxonomy_db`     | path | —            | Prebuilt MMseqs2 taxonomy database, required only when the above is `true`.                       |

## Output

| Path (relative to `--outdir`)                              | Content                                                                  |
| ---------------------------------------------------------- | ------------------------------------------------------------------------ |
| `assemblies/contigs/<sample>.contigs.fa.gz`                | Per-sample assembled contigs.                                            |
| `assemblies/bins/<sample>/<sample>_binN.fa`                | **Refined bins (MAGs) per sample**, from Binette.                        |
| `assemblies/bins/<sample>/quality_report.tsv`              | CheckM2 completeness and contamination per bin.                          |
| `assemblies/bins/<sample>/contig_to_bin.tsv`               | Which contig went into which bin.                                        |
| `binning/bams/`, `binning/bwa_index/`                      | Read mapping against each sample's own chromosomal contigs.              |
| `binning/metabat2/`, `binning/semibin2/`                   | The two upstream bin sets, before refinement.                            |
| `classification/per_contig/<sample>.tsv`                   | **Per-contig classification**, one row per contig (columns below).       |
| `classification/tiara/`                                    | Raw Tiara labels.                                                        |
| `classification/unbinned/`                                 | Contigs in no bin, and their taxonomy when enabled.                      |
| `classification/all.genes.classified.raw.tsv`              | Gene catalog with each gene's contig class attached.                     |
| `classification/all.genes.clusters.classified.refined.tsv` | **Gene clusters with a `multi_class` flag** — clusters spanning classes. |

The per-contig table carries `contig_id`, `sample`, `length`, `primary_class`, `classifier`, `lineage`, `confidence`, `bin_id`, `genomad_viral`, `genomad_plasmid` and `tiara_label`. `classifier` records which evidence source decided `primary_class`, so a disagreement between geNomad and Tiara is visible rather than silently resolved.

## Example

From reads:

```bash
nextflow run schirmer-lab/metagear-pipeline -profile docker \
  --workflow classification \
  --input clean.csv \
  --outdir results/ \
  --genomad_db /data/metagear/genomad \
  --checkv_db /data/metagear/checkv \
  --checkm2_db /data/metagear/checkm2
```

Reusing a prior `virus` run, so geNomad does not run twice:

```bash
nextflow run schirmer-lab/metagear-pipeline -profile docker \
  --workflow classification \
  --input clean.csv \
  --outdir results/ \
  --contigs_dir results/assemblies/contigs \
  --chromosome_dir results/virus/genomad/chromosome \
  --viral_ids_dir results/virus/per_sample \
  --plasmid_ids_dir results/virus/per_sample \
  --gene_clusters_tsv results/catalogs/genes/all.genes.clusters.tsv \
  --checkm2_db /data/metagear/checkm2
```

## Notes

- **Most contigs come back `unknown`, and that is expected.** A contig only gets a class if it lands in a bin or geNomad or Tiara calls it; short, low-coverage fragments do neither. Weight by length or coverage rather than counting rows before reading anything into the proportions.
- **`--enable_contig_taxonomy` is off deliberately.** Per-contig LCA on unbinned contigs is low-confidence, and against GTDB it needs upwards of 96 GiB and hours per sample. Bin-attributable contigs already get trustworthy lineage from `mag`; viral and plasmid contigs get theirs from geNomad. Turn it on for exploratory signal on the unbinned fraction, not for headline numbers.
- **The `biome` column selects SemiBin2's pretrained model.** It defaults to `global` when the column is absent or empty, and is passed through as SemiBin2's `--environment`.
- **Two binners, then refinement.** MetaBAT2 and SemiBin2 disagree often; Binette takes the union and keeps the best-scoring non-redundant set. A sample where one binner produces nothing still passes through.
- **`--chromosome_dir` does not check parameter consistency.** It trusts that the supplied partition came from a compatible geNomad run. If FDR thresholds or score calibration changed, re-run the detection rather than reusing.
- **Bins here are per-sample.** Dereplicating them into a cohort catalog with taxonomy is [mag](mag.md).
