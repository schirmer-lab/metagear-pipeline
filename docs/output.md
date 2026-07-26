# schirmer-lab/metagear-pipeline: Output

## Introduction

This document describes the output produced by the pipeline. The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## How the results tree is organised

MetaGEAR is a multi-workflow pipeline: a run produces only the directories belonging to the workflow you selected with `--workflow`. The layout is deliberately **additive** — every workflow writes into the same tree without overwriting another workflow's outputs, so you can point several runs at the same `--outdir` and accumulate results. This is what makes cross-workflow reuse (`--contigs_dir`, `--bins_dir`, or `--reuse-outputs` in the wrapper) work: a later workflow reads the artefacts an earlier one published.

Two conventions run through the tree. Directory names describe **contents** rather than the tool or operation that produced them, so `virus/detection/` rather than `virus/collect/`. Shareable files are named `<class>.<entity>.<role>` — for example `virus.contigs.representative.fa.gz` or `all.genes.rpkm.tsv` — so a file stays self-describing once it has been copied out of its directory. `<class>` is one of `all.genes`, `virus`, `plasmid`, `mag`, or `msp`.

Cohort-level deliverables cluster into three top-level directories shared by all analysis workflows:

| Directory     | Holds                                                                                             |
| ------------- | ------------------------------------------------------------------------------------------------- |
| `assemblies/` | Per-sample sequence products: contigs, assembly graphs, bins, dereplicated genomes                |
| `catalogs/`   | Cohort-level non-redundant reference sets, nested by granularity (genes, proteins, contigs, mags) |
| `abundance/`  | `<class>`-by-sample quantification matrices plus the BAMs and indexes used to produce them        |

## Quality control

<details markdown="1">
<summary>Output files</summary>

- `kneaddata/`
  - Host- and adapter-filtered read pairs, one set per sample.
  - Per-sample KneadData logs and read-count summaries.
- `multiqc/`
  - FastQC before and after trimming, aggregated across the cohort.

</details>

Produced by `qc_dna` and `qc_rna`. These reads are what every downstream analysis workflow expects as input. The `qc_rna` flow adds an rRNA depletion step but publishes to the same locations.

## Reference-based profiling

<details markdown="1">
<summary>Output files</summary>

- `metaphlan/`
  - `individual_profiles/<sample>_microbial_profile.txt`: per-sample MetaPhlAn 4 taxonomic profiles.
  - Merged cohort abundance table across all samples.
- `humann/`
  - Gene-family and pathway abundance and coverage tables, per-sample and merged.

</details>

Produced by `microbial_profiles`. This is the only analysis workflow that does not assemble, so it is substantially cheaper than the others and makes a good first pass on a new cohort.

## Assembly and gene catalog

<details markdown="1">
<summary>Output files</summary>

- `assemblies/contigs/<sample>.contigs.fa.gz`: per-sample MEGAHIT assemblies.
- `assemblies/assembly_graphs/<sample>.k119.fastg.gz`: assembly de Bruijn graphs, retained because third-party tools such as vRhyme and BinSPreader consume them.
- `prodigal/`: per-sample gene calls, raw and length-filtered.
- `catalogs/raw/`: pre-clustering concatenated gene and protein sets (intermediates).
- `catalogs/genes/`: `all.genes.representative.fa.gz` (the **representative gene catalog**) and `all.genes.clusters.tsv` (cluster membership).
- `catalogs/proteins/`: the equivalent protein-level catalog and cluster table.
- `clusters/`: cohort-wide gene-cluster annotations. These sit at top level rather than under `virus/` because they are not virus-specific; `classification` adds sibling files here without colliding.
- `abundance/all.genes/`: `bwa_index/`, `bams/<sample>.bam`, `per_batch/` provenance tables, and the merged `all.genes.{count,rpkm,tpm}.tsv` matrices.
- `annotations/amrfinder/`, `annotations/interproscan/`: functional annotation of the representative proteins.

</details>

Produced by `genes`. The three merged matrices under `abundance/all.genes/` are the headline deliverables — gene-by-sample counts, RPKM, and TPM. The `per_batch/` tables are kept as provenance for the merge and are not normally analysed directly.

## Viral and plasmid analysis

<details markdown="1">
<summary>Output files</summary>

- `virus/genomad/`, `virus/checkv/`: raw two-pass detection output (`virus/` and `provirus/` passes).
- `virus/per_sample/<sample>/`: per-sample detection intermediates with normalised names — `virus.summary.tsv`, `virus.ids.txt`, `virus.filtered.fna.gz`, plus plasmid counterparts.
- `virus/detection/`: cohort-level call catalogs, `<class>.calls.tsv` and `<class>.calls.filtered.tsv`.
- `virus/lifestyle/`: vOTU-level life-history annotations — AMG calls, iPHoP host predictions, and the AMG-to-catalog cross-reference.
- `virus/clusters/`: class-specific cluster artefacts, including `<class>-exclusive.genes.representative_ids.txt`.
- `virus/genes/<sample>/`, `virus/dramv/`, `virus/pharokka/`, `virus/virsorter2/`, `virus/iphop/`: supporting per-tool output.
- `catalogs/contigs/`: `virus.contigs.representative.fa.gz` (the **non-redundant vOTU catalog**) and the plasmid equivalent, each with its cluster TSV.
- `abundance/virus/`, `abundance/plasmid/`: class-specific abundance matrices in the same shape as `abundance/all.genes/`.

</details>

Produced by `virus`. This is the most expensive workflow in the pipeline; iPHoP host prediction dominates the runtime and is I/O bound rather than CPU bound.

## Contig classification and MAGs

<details markdown="1">
<summary>Output files</summary>

- `classification/per_contig/<sample>.contigs.tsv`: the **per-contig classification table**. Alongside the priority-driven `primary_class`, it keeps `genomad_viral`, `genomad_plasmid`, and `tiara_label` as separate evidence columns, so disagreements between callers stay visible rather than being collapsed.
- `classification/unbinned/`, `classification/tiara/`, `classification/contig_taxonomy/`: the unbinned partition, the eukaryote screen, and — only when `--enable_contig_taxonomy` is set — MMseqs2 per-contig LCA assignments.
- `assemblies/bins/<sample>/`: per-sample MAG FASTAs from Binette, plus `quality_report.tsv` and `contig_to_bin.tsv`.
- `binning/`: SemiBin2, MetaBAT2, and Binette process intermediates.
- `assemblies/dereplicated/`: cluster-winner genomes from dRep.
- `catalogs/mag/`: `mag.representative.fa.gz`, `contig_to_mag.tsv`, and `mag_catalog.csv` (per-cluster summary).
- `mag/drep/`, `mag/gtdbtk/`: dRep cluster tables and figures, and GTDB-Tk `bac120`/`ar53` taxonomy for the cluster representatives.
- `abundance/mag/`: MAG-by-sample abundance matrices.

</details>

Produced by `classification` and `mag`. `mag` builds on a completed `classification` run and reads its per-sample bins from disk. Per-contig taxonomy enrichment — joining the per-contig table against the dRep and GTDB-Tk outputs — is deliberately left out of the pipeline as a short pandas or DuckDB operation; see [the enrichment recipe](recipes/per-contig-taxonomy-enrichment.md).

## Metagenomic species pangenomes

<details markdown="1">
<summary>Output files</summary>

- `msp/mspminer/`: MSPminer co-abundance clustering output and MSP membership tables.
- `msp/gtdbtk/`: GTDB-Tk taxonomy assigned to MSP representative sequences.
- `abundance/msp/`: MSP-by-sample abundance matrices.

</details>

Produced by `msp`, which builds on a prior `genes` run. `msp` and `mag` are complementary cohort-level species catalogs — `msp` is gene-centric (co-abundance clustering over the gene catalog), `mag` is genome-centric (dRep over assembled bins). Each owns its own subtree, so running one never overwrites the other.

## Structural annotation

<details markdown="1">
<summary>Output files</summary>

- `annotations/phold/all.proteins.phold.tsv`: PHOLD structural-homology annotation per representative protein.
- `annotations/phold/virus.proteins.phold.tsv`: the viral subset, with a `relation` column recording whether each row was joined back to the viral catalog directly or propagated through the protein cluster table.
- `structures/phold_predict/`: persisted ProstT5 prediction shards.
- `structures/offsite_predict/`: only when `--structures_prepare_for_gpu` is set — a self-contained bundle (shards, runner script, README, checksum manifest) for running the GPU step on a separate machine.

</details>

Produced by `structures`, building on a prior `virus` or `genes` run. Predictions are persisted in the same on-disk shape whether they came from a local GPU run or an offsite one, so resuming with `--phold_predict_dir` behaves identically in both cases.

## MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory. The report also carries the consolidated software-versions table, so it doubles as a provenance record. For more information about how to use MultiQC reports, see <http://multiqc.info>.

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `<workflow>_mqc_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameters are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

The consolidated versions file is prefixed with the workflow name (`genes_mqc_versions.yml`, `virus_mqc_versions.yml`, and so on), so accumulating several workflows in one `--outdir` leaves one versions record per workflow rather than a single file each run overwrites. The `_mqc_` token is what makes MultiQC pick the file up as custom content and render it as the Software Versions section.

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
