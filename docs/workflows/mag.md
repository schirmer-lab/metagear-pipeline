# mag

Turn per-sample bins into one cohort MAG catalog: dereplicate across samples, place the representatives taxonomically, and quantify each MAG in every sample. Builds on a prior [classification](classification.md) run and reads its bins from disk — it does not bin anything itself.

## What it does

1. **Collect bins** — every `<sample>_binN.fa` under `--bins_dir` is gathered, along with the CheckM2 quality reports that `classification` published beside them.
2. **Pairwise ANI** — skani computes all-vs-all average nucleotide identity across the cohort's bins.
3. **Batching** — bins are partitioned into connected components at `--drep_batch_ani` (default 0.90) so dRep runs on tractable groups instead of one all-vs-all job.
4. **Dereplication** — dRep clusters at its usual 0.95 secondary threshold with skani as the backend, filtering on quality (`-comp 50 -con 10`), and picks one representative per cluster.
5. **Catalog assembly** — representatives are collected into a single FASTA with a contig-to-MAG map and a per-MAG summary.
6. **Taxonomy** — GTDB-Tk `classify_wf` places the representatives.
7. **Abundance** — reads from every sample are mapped to the MAG catalog with BWA and quantified with CoverM, giving a MAG-by-sample matrix.

## Inputs

- `--input` — the same samplesheet used for `classification`. Needed because abundance maps the original reads.
- `--bins_dir` — `assemblies/bins/` from the `classification` run. Auto-discovered by `--reuse-outputs`.
- `--gtdb_tk_db` — GTDB-Tk reference (R220 or compatible).

## Parameters

| Parameter           | Type  | Default      | Controls                                                                              |
| ------------------- | ----- | ------------ | --------------------------------------------------------------------------------------- |
| `--input`           | path  | _(required)_ | Samplesheet of clean reads.                                                           |
| `--outdir`          | path  | _(required)_ | Result directory.                                                                     |
| `--bins_dir`        | path  | _(required)_ | Per-sample Binette bins from `classification`.                                        |
| `--gtdb_tk_db`      | path  | —            | GTDB-Tk reference. Usually set once in `~/.metagear/metagear.config`.                 |
| `--drep_batch_ani`  | float | `0.90`       | ANI at which bins are grouped into dRep batches. See Notes before changing.           |

## Output

| Path (relative to `--outdir`)          | Content                                                          |
| -------------------------------------- | ------------------------------------------------------------------ |
| `catalogs/mag/mag.representative.fa.gz` | **The cohort MAG catalog** — one representative genome per cluster. |
| `catalogs/mag/mag_catalog.csv`         | **One row per MAG**: identity, source sample, quality.           |
| `catalogs/mag/contig_to_mag.tsv`       | Which contig belongs to which MAG.                               |
| `catalogs/mag/mag_contig_lengths.tsv`  | Contig lengths, for length-weighted statistics.                  |
| `mag/drep/Wdb.cohort.csv`              | dRep winners — the representative chosen per cluster.            |
| `mag/drep/Cdb.cohort.csv`              | Cluster membership: which bins collapsed together.               |
| `mag/drep/skani/`                      | Raw pairwise ANI.                                                |
| `mag/gtdbtk/cohort/`                   | **GTDB-Tk summaries** (`*.bac120.summary.tsv`, `*.ar53.summary.tsv`). |
| `abundance/mag/mag.count.tsv`          | **MAG-by-sample read counts.**                                   |
| `abundance/mag/mag.rpkm.tsv`           | **MAG-by-sample RPKM.**                                          |
| `abundance/mag/mag.tpm.tsv`            | **MAG-by-sample TPM.**                                           |

`Wdb.cohort.csv`, `mag_catalog.csv` and the number of representative FASTAs should all agree; a mismatch means dereplication and catalog assembly disagreed and is worth investigating before using the abundance matrices.

## Example

```bash
nextflow run schirmer-lab/metagear-pipeline -profile docker \
  --workflow mag \
  --input clean.csv \
  --outdir results/ \
  --bins_dir results/assemblies/bins \
  --gtdb_tk_db /data/metagear/gtdb_tk
```

With the wrapper, `--reuse-outputs` finds `--bins_dir` on its own:

```bash
metagear mag --input clean.csv --outdir results/ --reuse-outputs
```

## Notes

- **Batching does not change the answer.** Connected components at 0.90 ANI strictly contain dRep's 0.95 clustering threshold, so no pair that would have clustered together can be split across batches. The argument is written out in `bin/batch_bins.py`. Lowering `--drep_batch_ani` makes batches larger and slower; raising it past 0.95 breaks the guarantee.
- **Singletons still go through dRep.** A bin alone in its batch is not passed through untouched — it runs so the `-comp 50 -con 10` quality filter is applied uniformly across the cohort.
- **GTDB-Tk is the memory-heavy step.** It carries the `process_high_memory` label and is configured for 64 GB and 12 CPUs in `conf/resources.config`. Most of that is pplacer, which is largely single-threaded, so extra cores help the HMMER and ANI screening steps rather than the placement itself.
- **Abundance is competitive across the catalog.** Reads map against all MAGs at once, so closely related MAGs split reads between them rather than each claiming all. Read the matrices as relative within the cohort, not as absolute genome copy number.
- **Archaea are included.** GTDB-Tk places both domains; expect an `ar53` summary alongside `bac120`, often with very few rows.
