#!/usr/bin/env bash
# lib/common.sh

LAUNCH_DIR="$PWD"

declare -A require_input=(
    [download_databases]="false"
    [qc_dna]="true"
    [qc_rna]="true"
    [microbial_profiles]="true"
    [gene_call]="true"
)

# Prompt user with a default value if input is not provided.
function prompt_for_required_file() {
    local prompt_message="$1"

    read -p "$prompt_message: " input

    # Ask again if the input is empty or the file does not exist
    while [[ -z "$input" || ! -f "$input" ]]; do
        if [[ -f "$input" ]]; then
            break
        fi
        read -p "Invalid input. $prompt_message: " input
    done

    echo "${input}"
}


function run_workflows() {
    workflow="$1"
    shift

    local default_input_file="$LAUNCH_DIR/input_${workflow}.csv"
    local default_outdir="$LAUNCH_DIR/results"

    while [[ "$1" != "" ]]; do
        case "$1" in
            --input)
                shift
                input_file="$1"
                ;;
            --outdir)
                shift
                outdir="$1"
                ;;
            *)
                echo "Invalid option: $1"
                usage
                ;;
        esac
        shift
    done

    # Only execute if the workflow is in the require_input array
    if [[ "${require_input[$workflow]}" == "true" ]]; then
        if [ -z "$input_file" ]; then
            # Check if file.txt exists in the current directory
            if [ -f "$default_input_file" ]; then
                input_file="$default_input_file"
            else
                input_file=$(prompt_for_required_file "Please provide a valid input file")
                cp "$input_file" "$default_input_file"
            fi
        else
            cp "$input_file" "$default_input_file"
        fi
    fi

    if [ -z "$outdir" ]; then
        outdir="$default_outdir"
    fi

    mkdir -p "$outdir"

    echo "--workflow $workflow --input $default_input_file --outdir $outdir"
}

