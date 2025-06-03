#!/usr/bin/env bash

export NXF_SINGULARITY_CACHEDIR=/nfs/data/database/singularity

RUN_PROFILES="-profile singularity,singularity_custom"
NF_WORK="./nf_work"

nextflow run /nfs/home/emilio/.metagear/latest/main.nf \
        --workflow viral_analysis --input /nfs/arxiv/emilio/github/metagear-pipeline-dev/input_viral_analysis.csv --outdir /nfs/arxiv/emilio/github/metagear-pipeline-dev/results \
        -c /nfs/arxiv/emilio/github/metagear-pipeline-dev/.metagear/viral_analysis.config \
        $RUN_PROFILES -w \
        $NF_WORK -resume

