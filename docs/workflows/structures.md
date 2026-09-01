# structures

Recover function for proteins that sequence homology leaves unassigned, by comparing predicted structures instead of sequences. PHOLD translates each protein to the 3Di structural alphabet with ProstT5 and searches it against a curated database with Foldseek. Builds on a prior [genes](genes.md) or [virus](virus.md) run and reads the protein catalog from disk. Use it when a large fraction of your catalog is annotated only as hypothetical.

## What it does

1. **Scope selection** — the protein catalog is filtered to the subset worth annotating, controlled by `--structures_scope`. Viral proteins are topped up in every mode except `viral_only`, so viral coverage is guaranteed regardless of scope.
2. **Sharding** — the selected proteins are split into shards for the prediction step.
3. **Structure prediction** — ProstT5 converts each amino-acid sequence into a 3Di structural string. This is the expensive step and the one that benefits from a GPU.
4. **Merge** — per-shard predictions are combined into one 3Di FASTA.
5. **Structural search** — Foldseek compares the 3Di sequences against the PHOLD database, chunked so each invocation amortises its database load.
6. **Merge and split** — hits are merged into per-representative annotations, with the viral catalog emitted separately from the full one.

## Inputs

- `--input` — samplesheet. Used for provenance only; no reads are processed.
- `--representative_proteins`, `--representative_proteins_clusters`, `--viral_representative_proteins` — from the prior run. Discovered by `--reuse-outputs`.
- `--representative_proteins_annotations` — Pfam/InterProScan table. Optional: without it the scope filter treats every protein as unannotated.
- `--phold_db` — the PHOLD structural database (~7.7 GB), installed with `metagear download_databases --databases structures`.

## Parameters

| Parameter                               | Type | Default                | Controls                                                                 |
| --------------------------------------- | ---- | ---------------------- | ------------------------------------------------------------------------ |
| `--input`                               | path | _(required)_           | Samplesheet, for provenance.                                             |
| `--outdir`                              | path | _(required)_           | Result directory.                                                        |
| `--phold_db`                            | path | —                      | PHOLD structural database.                                               |
| `--representative_proteins`             | path | _(required)_           | Full protein catalog from `genes`/`virus`.                               |
| `--representative_proteins_clusters`    | path | _(required)_           | Protein cluster table, used to map annotations back to members.          |
| `--viral_representative_proteins`       | path | _(required)_           | Viral protein catalog, used for the guaranteed viral top-up.             |
| `--representative_proteins_annotations` | path | —                      | Pfam/InterProScan table. Without it every protein counts as unannotated. |
| `--structures_scope`                    | enum | `unannotated_plus_duf` | `all` \| `unannotated` \| `unannotated_plus_duf` \| `viral_only`.        |
| `--structures_use_cpu`                  | bool | `true`                 | CPU prediction. Set `false` on a GPU node for a 10–50× speedup.          |
| `--structures_batch_size`               | int  | `1`                    | ProstT5 batch size. Higher is faster on GPU and uses more VRAM.          |
| `--structures_predict_shard_size`       | int  | `5000`                 | Proteins per prediction shard.                                           |
| `--structures_compare_shard_size`       | int  | `100000`               | Proteins per Foldseek chunk.                                             |
| `--structures_prepare_for_gpu`          | bool | `false`                | Package the prediction step for an offsite GPU and stop. See Notes.      |
| `--phold_predict_dir`                   | path | —                      | Consume prediction output produced elsewhere, skipping prediction.       |

## Output

| Path (relative to `--outdir`)       | Content                                                                               |
| ----------------------------------- | ------------------------------------------------------------------------------------- |
| `structures/all.proteins.3di.fasta` | The 3Di structural translation of the selected proteins.                              |
| `structures/phold_predict/`         | Per-shard prediction output, and the resume point for offsite GPU runs.               |
| `structures/`                       | **Merged PHOLD annotations** for the full catalog and, separately, the viral catalog. |
| `structures/offsite_predict/`       | The self-contained bundle, only when `--structures_prepare_for_gpu` is set.           |

## Example

Default scope, on CPU:

```bash
nextflow run schirmer-lab/metagear-pipeline -profile docker \
  --workflow structures \
  --input clean.csv \
  --outdir results/ \
  --representative_proteins results/catalogs/proteins/all.proteins.representative.fa.gz \
  --representative_proteins_clusters results/catalogs/proteins/all.proteins.clusters.tsv \
  --viral_representative_proteins results/catalogs/proteins/virus.proteins.representative.fa.gz \
  --representative_proteins_annotations results/annotations/interproscan/all.proteins.FG_IPS_Pfam.tsv \
  --phold_db /data/metagear/phold
```

On a GPU node, larger batches:

```bash
metagear structures --input clean.csv --outdir results/ --reuse-outputs \
  --structures_use_cpu false --structures_batch_size 8
```

Offsite GPU, in two passes:

```bash
# 1. on the cluster — packages the work and stops
metagear structures --input clean.csv --outdir results/ --reuse-outputs \
  --structures_prepare_for_gpu true

# 2. rsync results/structures/offsite_predict/ to the GPU machine, run the bundled
#    script there, rsync the output back into results/structures/phold_predict/, then:
metagear structures --input clean.csv --outdir results/ --reuse-outputs \
  --phold_predict_dir results/structures/phold_predict
```

## Notes

- **Scope is the cost lever.** `all` annotates the whole catalog and is the most expensive thing in the pipeline. The default `unannotated_plus_duf` targets proteins with no Pfam hit plus those hitting only domains of unknown function — the proteins where a structural hit actually adds something. Annotating well-characterised proteins mostly reproduces what InterProScan already said.
- **Prediction dominates, and shard size decides the overhead.** ProstT5 reloads its model per shard at roughly 30–60 s each, so small shards spend their time loading. 5000 suits cohorts up to about a million proteins; raise it to 10000 at 500-sample scale.
- **Foldseek chunks work the opposite way.** Each `phold proteins-compare` pays a fixed database load, so chunks should be big enough to amortise it — aim for 30–60 minutes of work per chunk. Raise `--structures_compare_shard_size` for very large cohorts, lower it if you want finer-grained resume after a failure.
- **The offsite-GPU mode exists because the bundle needs no Nextflow.** The packaged script runs with plain Python and the model on a machine that may have neither containers nor a scheduler, which is what makes a borrowed GPU usable.
- **Viral coverage is guaranteed.** In every mode except `viral_only`, viral gene representatives not already in the chosen subset are added by a top-up step, so a scope chosen for cost reasons never silently drops the viral catalog.
- **Missing annotations degrade gracefully.** Without `--representative_proteins_annotations` the scope filter cannot tell annotated from unannotated and treats everything as unannotated, which is safe but larger than intended.
