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


# Run nextflow here... TODO
