# schirmer-lab/metagear-pipeline

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/schirmer-lab/metagear-pipeline)
[![GitHub Actions CI Status](https://github.com/schirmer-lab/metagear-pipeline/actions/workflows/nf-test.yml/badge.svg)](https://github.com/schirmer-lab/metagear-pipeline/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/schirmer-lab/metagear-pipeline/actions/workflows/linting.yml/badge.svg)](https://github.com/schirmer-lab/metagear-pipeline/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281%2Fzenodo.22233494-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.22233494)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**schirmer-lab/metagear-pipeline** is a bioinformatics pipeline for comprehensive metagenomic analysis. The pipeline processes shotgun metagenomic sequencing data through quality control, taxonomic profiling, functional annotation, and gene-centric analysis workflows.

> [!TIP]
> For easy installation, configuration, and usage, please refer to the **streamlined documentation** at **[MetaGEAR](https://metagear-platform.schirmerlab.de/)** and the **wrapper** at **[schirmer-lab/metagear-tools](https://github.com/schirmer-lab/metagear-tools)**.

The pipeline includes the following main workflows:

1. **Database Download** - Install required databases (Kneaddata, MetaPhlAn, HUMAnN)
2. **Quality Control** - DNA/RNA quality assessment and trimming
3. **Microbial Profiles** - Taxonomic and functional profiling using MetaPhlAn and HUMAnN
4. **Gene Analysis** - Gene-centric analysis workflow with optional contig catalogs

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

> [!IMPORTANT] > **For simplified usage and configuration**, please use the wrapper at **[schirmer-lab/metagear-tools](https://github.com/schirmer-lab/metagear-tools)**, documented at **[MetaGEAR](https://metagear-platform.schirmerlab.de/)**.

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

```bash
nextflow run schirmer-lab/metagear-pipeline \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

schirmer-lab/metagear-pipeline was written by Emilio Rios and Shen Jin.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

If you use schirmer-lab/metagear-pipeline for your analysis, please cite it using the following doi: [10.5281/zenodo.22233494](https://doi.org/10.5281/zenodo.22233494)

That DOI always resolves to the most recent release. To cite the exact version you ran, use the version-specific DOI shown on that record.

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file, grouped by the workflow that runs them. The MultiQC report produced by each run also contains a Methods Description section listing only the tools that particular run executed, with citations.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
