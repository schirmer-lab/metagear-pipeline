# qc_rna

Quality control for metatranscriptomic (RNA-seq) reads. Mechanically the same steps as [qc_dna](qc_dna.md), run against more reference databases: a host genome, usually a host transcriptome, and an rRNA reference.

## What it does

1. **FastQC** on the raw reads.
2. **Trim Galore** for adapter and quality trimming.
3. **KneadData** (v0.12.4) decontamination against every database in `--qc_rna_refdb`, in order.
4. **FastQC** again on the decontaminated reads.
5. **Read-count reporting** — per sample, per stage and per reference database.

## Inputs

- `--input` — CSV samplesheet (`sample,fastq_1,fastq_2`). One row per sample; see [usage.md](../usage.md).
- `--qc_rna_refdb` — directories holding Bowtie2 indexes. For metatranscriptomics, include at least a host genome **and** an rRNA reference; add the host transcriptome to remove residual host mRNA.

Without an rRNA reference the workflow still runs, but rRNA reads will dominate the output.

## Parameters

| Parameter               | Type  | Default      | Controls                                                                    |
| ----------------------- | ----- | ------------ | ----------------------------------------------------------------------------- |
| `--input`               | path  | _(required)_ | Samplesheet of raw FASTQ pairs.                                             |
| `--outdir`              | path  | _(required)_ | Where to write results.                                                     |
| `--qc_rna_refdb`        | array | `[]`         | Host genome, transcriptome and rRNA reference(s). Falls back to `kneaddata_refdb`. |
| `--qc_max_plot_samples` | int   | `100`        | Largest cohort for which the report draws the read-fate bar graph.          |

`qc_dna_refdb` and `qc_rna_refdb` are separate names so both assays can stay configured at once. Previously a single `kneaddata_refdb` in `metagear.config` had to be edited between a DNA and an RNA run, because the wrapper merges that file after the per-workflow config.

## Output

Same shape as `qc_dna`, with `qc_rna_` in place of `qc_dna_` on the report and the cohort table:

| Path (relative to `--outdir`)              | Content                                                         |
| ------------------------------------------ | ----------------------------------------------------------------- |
| `kneaddata/<sample>_paired_{1,2}.fastq.gz` | **Final clean reads** — host- and rRNA-depleted.                |
| `kneaddata/<sample>_kneaddata.log`         | Per-sample log, one decontamination block per reference.        |
| `kneaddata/stats/<sample>_read_counts.tsv` | Per-sample counts by stage, database and mate.                  |
| `kneaddata/stats/qc_rna_read_counts.tsv`   | The same for the whole cohort.                                  |
| `pipeline_info/qc_rna_multiqc_report.html` | Consolidated MultiQC report.                                    |

Expect the final read count to drop substantially compared with DNA. That is the rRNA depletion working.

## Why the per-database breakdown matters here

For DNA there is usually one reference, so "contamination" is unambiguous. For RNA there are several, and a single combined figure cannot tell you whether a sample lost its reads to host genomic DNA, host mRNA, or rRNA — which are three different problems with three different fixes.

The report therefore carries one series per reference database rather than one "contamination" bucket. Reading it:

- a large **rRNA** fraction is expected, and is the depletion protocol's signature;
- a large **host transcriptome** fraction points at host RNA carry-over;
- a large **host genome** fraction points at DNA contamination of the RNA prep.

See [qc_dna](qc_dna.md#serial-and-parallel-filtering) for how the per-database counts are derived and why the pipeline pins serial filtering. It matters more for RNA than for DNA: with several references a read can match more than one, and only serial filtering attributes it to exactly one.

## Example

```bash
nextflow run schirmer-lab/metagear-pipeline -profile singularity \
  --workflow qc_rna \
  --input raw_rna_samples.csv \
  --outdir qc_rna/ \
  -c my.config   # sets params.qc_rna_refdb = ["/data/hg38_bt2", "/data/hg38_mRNA_bt2", "/data/SILVA_128_rRNA_bt2"]
```

## Notes

- **Expect heavy attrition.** A typical RNA library can be 70–90% rRNA; the cleaned output may be a fraction of the raw input. This is by design, and the per-database breakdown is how you confirm it went where you expected.
- **Database order is the filtering order.** Under serial filtering each reference sees what the previous one left, so the counts depend on the order you list them in. Keep the order stable across a cohort or the columns are not comparable between samples.
- **Downstream compatibility.** The output FASTQs are interchangeable with `qc_dna` output. Remember that abundance from RNA reflects activity, not presence.
