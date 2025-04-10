#!/bin/bash

# Default values
module=""
subworkflow=""
workflow=""
config=""

# Support vairables
script_dir="$(cd "$(dirname "$0")" && pwd)"
launch_dir="$PWD"

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
cat <<EOF > $launch_dir/metagear_entrypoint.config
params {
    module = "$module"
    subworkflow = "$subworkflow"
    workflow = "$workflow"
}
EOF

# Create user config file if it does not exist
if [ ! -f "$launch_dir/metagear_user.config" ]; then
    cp $script_dir/metagear.config $launch_dir/metagear_user.config
    echo ""
    echo "Configuration file was not found. New file is created.."
    echo "Please edit this file before continuing: $launch_dir/metagear_user.config"
    echo ""
    exit 0
fi


if [ ! -f "$launch_dir/metagear_run.config" ]; then
    additional_files=(  $script_dir/conf/metagear/$config $script_dir/metagear.config $launch_dir/metagear_user.config $launch_dir/metagear_entrypoint.config )
    config_files=( $script_dir/conf/metagear/*.config )
    all_files=( "${config_files[@]}" "${additional_files[@]}" )
    $script_dir/metagear_configure.sh ${all_files[@]}
fi



more_args=""
# Dummy input if workflow is setup
if [[ "$workflow" == "setup" ]]; then
    temp_file=$(mktemp)
    echo "sample,fastq_1,fastq_2" > $temp_file
    more_args="--input=$temp_file"
fi

if [ ! -f "$launch_dir/metagear_run.sh" ]; then
    $script_dir/metagear_initialize.sh $launch_dir $script_dir
fi

echo "Starting Nextflow pipeline..."
$launch_dir/metagear_run.sh $more_args -resume

