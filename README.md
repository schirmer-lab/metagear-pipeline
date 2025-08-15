# schirmer-lab/metagear

[![GitHub Actions CI Status](https://github.com/schirmer-lab/metagear-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/schirmer-lab/metagear-pipeline/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/schirmer-lab/metagear-pipeline/actions/workflows/linting.yml/badge.svg)](https://github.com/schirmer-lab/metagear-pipeline/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**schirmer-lab/metagear** is a bioinformatics pipeline for comprehensive metagenomic analysis. The pipeline processes shotgun metagenomic sequencing data through quality control, taxonomic profiling, functional annotation, and gene-centric analysis workflows.

> [!TIP]
> For easy installation, configuration, and usage, please refer to the **streamlined documentation and wrapper** at: **[schirmer-lab/metagear](https://schirmer-lab.github.io/metagear)**

The pipeline includes the following main workflows:

1. **Database Download** - Install required databases (Kneaddata, MetaPhlAn, HUMAnN)
2. **Quality Control** - DNA/RNA quality assessment and trimming
3. **Microbial Profiles** - Taxonomic and functional profiling using MetaPhlAn and HUMAnN
4. **Gene Analysis** - Gene-centric analysis workflow with optional contig catalogs

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

> [!IMPORTANT] > **For simplified usage and configuration**, please use the wrapper and documentation at: **[schirmer-lab/metagear](https://github.com/schirmer-lab/metagear)**

### Quick Start

Prepare a samplesheet with your input data:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
SAMPLE1,sample1_R1.fastq.gz,sample1_R2.fastq.gz
SAMPLE2,sample2_R1.fastq.gz,sample2_R2.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run schirmer-lab/metagear \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

schirmer-lab/metagear was originally written by Shen Jin, Emilio Rios, Svenja Schorlemmer.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use schirmer-lab/metagear for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
