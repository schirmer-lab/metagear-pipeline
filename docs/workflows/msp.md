# msp

Group the gene catalog into MetaSpecies Pangenomes: sets of genes that rise and fall together across the cohort, which in practice correspond to species-level units. Builds on a prior [genes](genes.md) run and reads its catalog and abundance matrices from disk. Use it when you want species-level structure derived from your own data rather than from a reference.

## What it does

1. **Co-abundance clustering** — MSPminer groups genes whose abundance covaries across samples into MSPs, separating each into a core set and accessory sets.
2. **Pangenome sequences** — the genes of each MSP are collected into a per-MSP FASTA.
3. **MSP abundance** — each MSP gets a per-sample abundance from the median of its core genes, which is robust to individual genes mapping badly.
4. **Taxonomy** — GTDB-Tk places the MSP representative sequences.
5. **MetaPhlAn cross-walk** — MSPs are regressed against MetaPhlAn profiles so each MSP can be related to a named reference taxon where one exists.

## Inputs

- `--input` — the same samplesheet as the `genes` run. Needed so MetaPhlAn profiles can be computed if they were not supplied.
- `--representative_genes`, `--representative_genes_count`, `--representative_genes_rpkm` — from the prior `genes` run. All three are required; `--reuse-outputs` discovers them.
- `--gtdb_tk_db` — GTDB-Tk reference.

## Parameters

| Parameter                     | Type | Default      | Controls                                                                       |
| ----------------------------- | ---- | ------------ | -------------------------------------------------------------------------------- |
| `--input`                     | path | _(required)_ | Samplesheet of clean reads.                                                    |
| `--outdir`                    | path | _(required)_ | Result directory.                                                              |
| `--representative_genes`      | path | _(required)_ | Representative gene catalog from `genes`.                                      |
| `--representative_genes_count`| path | _(required)_ | Gene-by-sample count matrix.                                                   |
| `--representative_genes_rpkm` | path | _(required)_ | Gene-by-sample RPKM matrix.                                                    |
| `--metaphlan_profiles`        | path | —            | Pre-computed MetaPhlAn profiles. When omitted they are computed from the reads.|
| `--gtdb_tk_db`                | path | —            | GTDB-Tk reference. Usually set once in `~/.metagear/metagear.config`.          |

The three catalog inputs are marked optional in the schema because `--reuse-outputs` supplies them, but the workflow exits early with a clear message if they are neither passed nor discovered.

## Output

| Path (relative to `--outdir`)          | Content                                                                   |
| -------------------------------------- | --------------------------------------------------------------------------- |
| `msp/mspminer/mspminer/all_msps.tsv`   | **MSP membership** — every gene, its MSP, and whether it is core or accessory. |
| `msp/mspminer/mspminer/all_core_seeds.tsv` | The core gene set that defines each MSP.                                |
| `msp/mspminer/mspminer/genes_filtered_out.tsv` | Genes MSPminer excluded, with the reason.                            |
| `msp/pangenome_sequences/`             | One FASTA per MSP.                                                        |
| `msp/gtdbtk/pangenome/`                | **GTDB-Tk taxonomy** for the MSP representatives.                         |
| `msp/msp_metaphlan4_LM.bestR2.txt`     | **Best MetaPhlAn match per MSP**, by regression R².                       |
| `msp/msp_metaphlan4_LM.full.txt`       | All MSP-to-taxon regressions, not just the best.                          |
| `abundance/msp/msp_abundance.median.tsv` | **MSP-by-sample abundance**, median over core genes.                    |

## Example

```bash
nextflow run schirmer-lab/metagear-pipeline -profile docker \
  --workflow msp \
  --input clean.csv \
  --outdir results/ \
  --representative_genes results/catalogs/genes/all.genes.representative.fa.gz \
  --representative_genes_count results/abundance/all.genes/all.genes.count.tsv \
  --representative_genes_rpkm results/abundance/all.genes/all.genes.rpkm.tsv \
  --gtdb_tk_db /data/metagear/gtdb_tk
```

With the wrapper:

```bash
metagear msp --input clean.csv --outdir results/ --reuse-outputs
```

## Notes

- **Co-abundance needs samples, not depth.** MSPminer infers structure from how genes covary across the cohort, so the number of samples is what it has to work with — sequencing any one sample more deeply does not add signal here.
- **MSPs are not species, and the cross-walk is a hint.** An MSP is a set of co-varying genes; the MetaPhlAn linkage is a regression, not an assignment. Check the R² in `bestR2` before naming an MSP after a taxon.
- **The catalog must match the abundance matrices.** MSPminer keys on gene identifiers, so the catalog and the count/RPKM tables have to come from the same `genes` run. If the catalog was rebuilt in between — re-clustering picks different representatives — the identifiers no longer line up and most genes drop out. `--reuse-outputs` takes all three from one output tree, which is the safe path.
- **Core versus accessory is the useful distinction.** Core genes are present across the samples where the MSP is present and are what the abundance is computed from; accessory genes vary between strains and are where the interesting biology usually is.
- **This is the gene-centric counterpart to [mag](mag.md).** MSPs come from abundance covariation across the cohort; MAGs come from assembly and binning within each sample. They answer similar questions by independent routes and are worth comparing.
