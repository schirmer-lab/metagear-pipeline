# download_databases

One-time setup of the reference databases consumed by every other workflow. Run this once per installation (and re-run when you want to refresh a database).

## What it does

Downloads and prepares the following references, sized at tens of gigabytes each:

**For `microbial_profiles` and `genes`:**

- **KneadData host reference** (e.g. human genome) for host-read decontamination during QC.
- **MetaPhlAn 4 marker database** for taxonomic profiling.
- **HUMAnN 3 ChocoPhlAn** (nucleotide) and **UniRef90** (protein) databases for functional profiling.
- **GTDB-Tk** reference for prokaryotic taxonomic placement of MSPs (gene-analysis only).

**For `virus`:**

- **geNomad**, **CheckV**, **VirSorter2** databases for viral/plasmid detection.
- **Pharokka** database for prokaryotic-virus gene annotation.
- **DRAM-v** database for auxiliary metabolic gene (AMG) discovery.
- **iPHoP** database for host prediction of viral contigs.
- **AMRFinderPlus** database for antimicrobial-resistance gene annotation.

Each tool is downloaded by its native installer wrapped in a Nextflow module; nothing is repackaged.

## Inputs

This workflow does not read a samplesheet — there is no `--input` argument. Instead, you tell it where to put each database via the path parameters below. The `metagear-tools` CLI manages these paths for you via its config; if you run Nextflow directly you must pass them yourself (typically through a `-params-file`).

## Parameters

| Parameter            | Type   | Default      | Controls                                                       |
| -------------------- | ------ | ------------ | -------------------------------------------------------------- |
| `--databases`        | string | `all`        | Which database set to download: `all`, `genes`, or `virus`.    |
| `--outdir`           | path   | _(required)_ | Output directory; database files land in `<outdir>/<tool>/...` |
| `--kneaddata_refdb`  | array  | `[""]`       | Destination(s) for KneadData host reference(s).                |
| `--metaphlan_db`     | path   | —            | Destination for MetaPhlAn 4 database.                          |
| `--humann3_nucleo`   | path   | —            | Destination for HUMAnN 3 ChocoPhlAn nucleotide database.       |
| `--humann3_uniref90` | path   | —            | Destination for HUMAnN 3 UniRef90 protein database.            |
| `--gtdb_tk_db`       | path   | —            | Destination for GTDB-Tk reference.                             |
| `--genomad_db`       | path   | —            | Destination for geNomad database.                              |
| `--checkv_db`        | path   | —            | Destination for CheckV database.                               |
| `--virsorter2_db`    | path   | —            | Destination for VirSorter2 database.                           |
| `--pharokka_db`      | path   | —            | Destination for Pharokka database.                             |
| `--dram_db`          | path   | —            | Destination for DRAM-v database.                               |
| `--iphop_db`         | path   | —            | Destination for iPHoP database.                                |
| `--amrfinder_db`     | path   | —            | Destination for AMRFinderPlus database.                        |

## Output

Each tool writes into the directory you supplied for its `*_db` parameter. The pipeline does not republish those files into `<outdir>` (the modules' `publishDir` is intentionally suppressed) — the databases stay where the installers placed them.

Typical layout (paths are whatever you passed):

```
<metaphlan_db>/                      # MetaPhlAn 4 marker DB
<humann3_nucleo>/                    # ChocoPhlAn pangenomes
<humann3_uniref90>/                  # UniRef90 protein DB
<kneaddata_refdb[0]>/                # Bowtie2 host genome index
<gtdb_tk_db>/                        # GTDB-Tk release
<genomad_db>/                        # geNomad models
<checkv_db>/                         # CheckV references
<pharokka_db>/                       # Pharokka databases
<dram_db>/                           # DRAM-v references
<iphop_db>/                          # iPHoP host genomes
<amrfinder_db>/                      # AMRFinderPlus DB

<outdir>/pipeline_info/              # Nextflow execution metadata
```

## Example

```bash
# Download everything (recommended on first install)
nextflow run schirmer-lab/metagear-pipeline -profile docker \
  --workflow download_databases \
  --databases all \
  --outdir databases/ \
  -params-file db-paths.yaml
```

`db-paths.yaml`:

```yaml
metaphlan_db: /data/metagear/metaphlan
humann3_nucleo: /data/metagear/humann/chocophlan
humann3_uniref90: /data/metagear/humann/uniref90_diamond
kneaddata_refdb: [/data/metagear/kneaddata/Homo_sapiens]
gtdb_tk_db: /data/metagear/gtdb_tk
genomad_db: /data/metagear/genomad
checkv_db: /data/metagear/checkv
pharokka_db: /data/metagear/pharokka
dram_db: /data/metagear/dram
iphop_db: /data/metagear/iphop
amrfinder_db: /data/metagear/amrfinder
```

## Notes

- **Disk and time.** A full download is in the order of 100+ GB and several hours, mostly bound by external download speed.
- **Run once, point everything at it.** Downstream workflows read from these same paths via the same parameters, so once databases are in place all subsequent runs only need `--input`, `--outdir`, and `-profile`.
- **Selective re-runs.** Use `--databases genes` or `--databases virus` to refresh only one set without re-downloading the other.
- **HUMAnN config.** The HUMAnN installer is invoked with `--update-config no`; HUMAnN's global config is not modified on the host.
- **Wrapper convenience.** [`metagear-tools`](https://github.com/schirmer-lab/metagear-tools) ships templates that wire all of these paths into a single user config — recommended over hand-writing the params file above.
