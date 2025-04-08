#!/bin/bash

# Default values
module=""
subworkflow=""
workflow=""
config=""

# Parse the single command-line argument.
# Only one of these is expected.
for arg in "$@"; do
    case $arg in
        --module=*)
            module="${arg#*=}"
            config="${module}.config"
            ;;
        --subworkflow=*)
            subworkflow="${arg#*=}"
            config="${subworkflow}.config"
            ;;
        --workflow=*)
            workflow="${arg#*=}"
            config="${workflow}.config"
            ;;
        *)
            echo "Usage: $0 [--module=value | --subworkflow=value | --workflow=value]"
            exit 1
            ;;
    esac
done

# Write the structure with the provided value into params.txt
cat <<EOF > metagear_run.config
params {
    module = "$module"
    subworkflow = "$subworkflow"
    workflow = "$workflow"
}
EOF

additional_files=(  conf/metagear/$config metagear_user.config metagear_run.config )
config_files=( conf/metagear/*.config )
all_files=( "${config_files[@]}" "${additional_files[@]}" )

./metagear_configure.sh ${all_files[@]}



# export REPO=/nfs/arxiv/emilio/github/metagear-pipeline
# export NXF_SINGULARITY_CACHEDIR=/nfs/data/database/singularity

# INPUT_FILE=/nfs/data/projects/Core_Pipelines/configurations/input/kch_pilot/clean_for_metaphlan.csv
# RESULTS_DIR=/nfs/arxiv/emilio/results/metagear
# WORK_DIR=/nfs/arxiv/emilio/nf_work/metagear

# mkdir -p $RESULTS_DIR $WORK_DIR && cd $WORK_DIR

# clear && nextflow $REPO/main.nf -c $REPO/metagear_run.config -profile singularity,singularity_schirmer -resume -w $WORK_DIR
