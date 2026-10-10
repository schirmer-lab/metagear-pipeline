# qc_dna

Quality control for shotgun-metagenomic DNA reads. Produces decontaminated, adapter-trimmed paired reads that are ready for any of the analysis workflows.

## What it does

1. **FastQC** on the raw reads.
2. **Trim Galore** for adapter and quality trimming. Process arg: `--phred33 --quality 0 --stringency 5 --length 10` (light defaults; quality trimming is delegated to KneadData).
3. **KneadData** (v0.12.4, internally running Trimmomatic and Bowtie2) for host decontamination against every database in `--qc_dna_refdb`.
4. **FastQC** again on the decontaminated reads.
5. **Read-count reporting** — per sample, per stage and per reference database.

`qc_rna` runs these same steps against a different set of references. See [qc_rna](qc_rna.md).

## Inputs

- `--input` — CSV samplesheet (`sample,fastq_1,fastq_2`). One row per sample; see [usage.md](../usage.md) for the contract.
- `--qc_dna_refdb` — one or more directories, each holding a Bowtie2 index to decontaminate against. For DNA this is normally a single host genome.

Mate files may be named `_1`/`_2`, `_R1`/`_R2` or `_R1_001`/`_R2_001`. Which file is which mate comes from the samplesheet columns, not the file name.

## Parameters

| Parameter               | Type  | Default      | Controls                                                                                            |
| ----------------------- | ----- | ------------ | --------------------------------------------------------------------------------------------------- |
| `--input`               | path  | _(required)_ | Samplesheet of raw FASTQ pairs.                                                                     |
| `--outdir`              | path  | _(required)_ | Where to write results.                                                                             |
| `--qc_dna_refdb`        | array | `[]`         | Reference database(s) for decontamination. Falls back to `kneaddata_refdb`.                         |
| `--qc_max_plot_samples` | int   | `100`        | Largest cohort for which the report draws the read-fate bar graph. Beyond it, the table is written. |

`qc_dna_refdb` and `qc_rna_refdb` are separate names on purpose. The `metagear` wrapper merges your `metagear.config` *after* the per-workflow config, so a single shared `kneaddata_refdb` there would override whichever assay you selected. With one name per assay, both stay configured and neither has to be edited between runs.

## Output

| Path (relative to `--outdir`)              | Content                                                     |
| ------------------------------------------ | ------------------------------------------------------------- |
| `kneaddata/<sample>_paired_{1,2}.fastq.gz` | **Final clean reads.** The input to every other workflow.   |
| `kneaddata/<sample>_kneaddata.log`         | Per-sample KneadData log.                                   |
| `kneaddata/stats/<sample>_read_counts.tsv` | Per-sample counts by stage, database and mate.              |
| `kneaddata/stats/qc_dna_read_counts.tsv`   | The same for the whole cohort, one row per count.           |
| `fastqc/`                                  | FastQC before and after cleaning.                           |
| `pipeline_info/qc_dna_multiqc_report.html` | Consolidated MultiQC report.                                |

The clean reads are the only output to keep long-term; the rest is diagnostic.

## Reading the read counts

KneadData counts reads at every stage and once per reference database, and writes those counts to its log. The pipeline parses them into `*_read_counts.tsv` and renders them in MultiQC as:

- a **stacked bar graph**, one bar per sample, splitting raw reads into those lost to trimming, those matched by each database, and those surviving;
- a **table** carrying the same numbers, sortable and exportable from the report.

### Serial and parallel filtering

How a per-database count should be read depends on how KneadData filtered:

- **Serial** (what this pipeline pins, via `--serial`). Each database filters what the previous one left. A read matching two databases is attributed to the first, removals **partition** the reads, and the parts sum to the raw count. This is what makes the stacked bar graph meaningful.
- **Parallel.** Every database filters the trimmed reads independently and the results are intersected. A read matching two databases is counted against both, so removals **overlap** and sum to more than the whole.

The parser reads the mode from the log rather than assuming it, because KneadData 0.12 forces serial whenever `--processes` is 1. Under parallel filtering the bar graph is suppressed rather than drawn misleadingly; the table still carries every per-database column.

### Large cohorts

A bar per sample stops being readable somewhere around a hundred samples. Past `--qc_max_plot_samples` the bar graph is dropped and only the table is written. The table is sortable and exportable but not paginated; above MultiQC's own row limit it renders as a violin plot instead of rows.

`kneaddata/stats/qc_dna_read_counts.tsv` is written on every run regardless, so the report's choice of plot never limits what can be analysed.

## Example

```bash
nextflow run schirmer-lab/metagear-pipeline -profile singularity \
  --workflow qc_dna \
  --input raw_samples.csv \
  --outdir qc/ \
  -c my.config   # sets params.qc_dna_refdb = ["/data/metagear/hg38_bt2"]
```

After this, build a `clean.csv` pointing at `kneaddata/<sample>_paired_{1,2}.fastq.gz` and feed it to `microbial_profiles`, `genes`, or `virus`.

## Notes

- **Trim Galore is intentionally light** (`--quality 0`); the heavier quality trimming is done by KneadData's internal Trimmomatic step with `SLIDINGWINDOW:4:15 MINLEN:50`. Don't tighten Trim Galore unless you know why.
- **`HEADCROP:15` removes the first 15 bases** of each read inside KneadData. Appropriate for typical Illumina libraries with biased initial bases; drop it via a custom config if your library doesn't have that artefact.
- **Single-end input** is accepted — leave `fastq_2` empty in the samplesheet.
- **FASTQ headers are not rewritten by the pipeline.** KneadData 0.12 normalises read identifiers itself, truncating at the first space and appending `#0/1` or `#0/2` by mate. The pipeline's own header-rewriting step existed for KneadData 0.10, which did not, and has been removed.
