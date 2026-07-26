# microbial_profiles

Reference-based taxonomic and functional profiling of microbial communities. Use this when you want to know "which species are here and what can they do" without committing to a de novo assembly.

## What it does

1. **MetaPhlAn 4** (v4.1.1) maps reads against species-specific marker genes and produces a relative-abundance profile per sample. With `--include_virus true` the `--profile_vsc` flag is added and a parallel viral profile is emitted.
2. **MetaPhlAn merge** combines per-sample profiles into a single abundance table across all samples.
3. **HUMAnN 3** uses the MetaPhlAn output to pick species-specific pangenomes from ChocoPhlAn, aligns reads against them, and falls back to UniRef90 for unmapped reads. Output is per-sample gene-family and pathway-abundance tables (normalized to CPM via `--units cpm`).
4. **HUMAnN merge** combines per-sample gene-family and pathway tables.

You can skip the MetaPhlAn step entirely with `--metaphlan_dir <dir>` if you already have profiles from a prior run.

## Inputs

- `--input` — samplesheet of QC'd reads.
- `--metaphlan_db` — MetaPhlAn 4 marker database (from [download_databases](download_databases.md)).
- `--humann3_nucleo` — HUMAnN ChocoPhlAn nucleotide database.
- `--humann3_uniref90` — HUMAnN UniRef90 protein database.
- _Optional:_ `--metaphlan_dir <dir>` — directory of pre-computed `*_microbial_profile.txt` files, one per sample; the MetaPhlAn step is skipped.

## Parameters

| Parameter            | Type    | Default      | Controls                                               |
| -------------------- | ------- | ------------ | ------------------------------------------------------ |
| `--input`            | path    | _(required)_ | Samplesheet of clean reads.                            |
| `--outdir`           | path    | _(required)_ | Result directory.                                      |
| `--metaphlan_db`     | path    | —            | MetaPhlAn 4 database.                                  |
| `--humann3_nucleo`   | path    | —            | ChocoPhlAn DB.                                         |
| `--humann3_uniref90` | path    | —            | UniRef90 DB.                                           |
| `--metaphlan_dir`    | path    | _(unset)_    | Skip MetaPhlAn and reuse profiles from this directory. |
| `--include_virus`    | boolean | `false`      | Also produce viral marker profiles (`--profile_vsc`).  |

## Output

| Path (relative to `--outdir`)                                  | Content                                                          |
| -------------------------------------------------------------- | ---------------------------------------------------------------- |
| `metaphlan/individual_profiles/<sample>_microbial_profile.txt` | Per-sample MetaPhlAn species table (relative abundance).         |
| `metaphlan/merged_microbial_profiles.txt`                      | **Species-by-sample abundance matrix** across the whole study.   |
| `metaphlan/merged_viral_profiles.txt`                          | Same shape, viral markers. Only present with `--include_virus`.  |
| `humann/<sample>_genefamilies.tsv`                             | Per-sample gene-family abundance (UniRef90 IDs, CPM-normalized). |
| `humann/<sample>_pathabundance.tsv`                            | Per-sample MetaCyc pathway abundance (CPM-normalized).           |
| `humann/merged_genefamilies.tsv`                               | **Gene-family-by-sample abundance matrix.**                      |
| `humann/merged_pathabundance.tsv`                              | **Pathway-by-sample abundance matrix.**                          |
| `pipeline_info/microbial_profiles_multiqc_report.html`         | Consolidated MultiQC.                                            |

The two `merged_*` tables in `humann/` and the `merged_microbial_profiles.txt` in `metaphlan/` are the typical inputs to downstream statistical analysis.

## Example

```bash
nextflow run schirmer-lab/metagear-pipeline -profile docker \
  --workflow microbial_profiles \
  --input clean.csv \
  --outdir profiles/ \
  --metaphlan_db /data/metagear/metaphlan \
  --humann3_nucleo /data/metagear/humann/chocophlan \
  --humann3_uniref90 /data/metagear/humann/uniref90_diamond
```

Reusing existing MetaPhlAn profiles (e.g. from a prior `genes` run that emitted them):

```bash
nextflow run schirmer-lab/metagear-pipeline -profile docker \
  --workflow microbial_profiles \
  --input clean.csv \
  --outdir profiles/ \
  --metaphlan_dir prior_run/metaphlan/individual_profiles/ \
  --humann3_nucleo /data/metagear/humann/chocophlan \
  --humann3_uniref90 /data/metagear/humann/uniref90_diamond
```

## Notes

- **MetaPhlAn 4 is reference-based.** Species that are not in the marker database will not be reported, even if they are abundant in your sample. For novel-organism discovery, complement with [genes](genes.md).
- **HUMAnN memory.** HUMAnN scales with diversity; 16 GB+ per sample is typical and large/diverse samples (e.g. soil) may need much more.
- **CPM units.** `humann_renorm_table` is run with `--units cpm`. If you need RPK or relative abundance, post-process the merged tables with the HUMAnN utilities.
- **Viral profiling is optional and reference-limited.** `--include_virus true` produces MetaPhlAn's VSC table; for de novo viral discovery use [virus](virus.md) instead.
