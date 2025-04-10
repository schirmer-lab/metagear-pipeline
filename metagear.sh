#!/bin/bash

# Default values
module=""
subworkflow=""
workflow=""
config=""

# Support vairables
script_dir="$(cd "$(dirname "$0")" && pwd)"
calling_dir="$PWD"

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
cat <<EOF > $calling_dir/metagear_run.config
params {
    module = "$module"
    subworkflow = "$subworkflow"
    workflow = "$workflow"
}
EOF

additional_files=(  $script_dir/conf/metagear/$config $script_dir/metagear.config $calling_dir/metagear_user.config $calling_dir/metagear_run.config )
config_files=( $script_dir/conf/metagear/*.config )
all_files=( "${config_files[@]}" "${additional_files[@]}" )

./metagear_configure.sh ${all_files[@]}

$more_args=""
# Dummy input if workflow is setup
if [[ "$workflow" == "setup" ]]; then
    temp_file=$(mktemp)
    echo "sample,fastq_1,fastq_2" > $temp_file
    more_args="--input=$temp_file"
fi


nextflow run $script_dir/main.nf -c $calling_dir/metagear_run.config -profile docker $more_args
