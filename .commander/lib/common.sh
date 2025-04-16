#!/usr/bin/env bash
# lib/common.sh

SCRIPT_DIR="$(dirname "$(realpath "$0")")"

source $SCRIPT_DIR/lib/system_utils.sh


declare -A commands=(
    [download_databases]="Install Databases (Kneaddata, Metaphlan, Humann)"
    [qc_dna]="Quality Control for DNA"
    [qc_rna]="Quality Control for RNA"
    [microbial_profiles]="Get microbial profiles with Metaphlan and Humann"
)


# Usage message
function usage() {
    echo ""
    echo "Usage: metagear <command> [options]"
    echo "Commands:"
    for cmd in "${!commands[@]}"; do
        printf "  %-20s %s\n" "$cmd" "${commands[$cmd]}."
    done
    echo ""
    exit 1
}


function check_command {
    # Check if the command exists in the commands array
    if [[ -z "${commands[$1]}" ]]; then
        echo "Error: Command '$1' not found."
        usage
        exit 1
    fi
}


check_requirements() {
    # Array to store error messages.
    local errors=()

    # Check Bash version 4+.
    if [ "${BASH_VERSINFO[0]}" -lt 4 ]; then
        errors+=("Bash version 4 or higher is required (found version ${BASH_VERSINFO[0]}).")
    fi

    # Check for nextflow.
    if ! command -v nextflow >/dev/null 2>&1; then
        errors+=("Nextflow is not installed.")
    fi

    # Check for a container engine: either singularity or docker.
    if ! command -v singularity >/dev/null 2>&1 && ! command -v docker >/dev/null 2>&1; then
        errors+=("Neither Singularity nor Docker is installed (one is required).")
    fi

    # If there are missing requirements, report them and exit.
    if [ ${#errors[@]} -gt 0 ]; then
        echo "The following requirements are missing:" >&2
        for error in "${errors[@]}"; do
            echo " - $error" >&2
        done
    fi

}


function check_metagear_home() {

    user_config_file=$HOME/.metagear/metagear.config
    user_env_file=$HOME/.metagear/metagear.env

    if [ ! -d "$HOME/.metagear" ]; then

        mkdir -p "$HOME/.metagear"

        echo "System resources:"
        echo "-----------------"

        total_cpu_count=$(get_cpu_count)
        echo "CPU Count: ${total_cpu_count}"

        total_memory_gb=$(get_total_memory_gb)
        echo "Installed RAM: ${total_memory_gb} GB"

        cp $1/templates/metagear.config $user_config_file

        sed -i "s/max_memory = '[^']*GB'/max_memory = '${total_memory_gb}GB'/" "$user_config_file"
        sed -i "s/max_cpus = [^']*/max_cpus = ${total_cpu_count}/" "$user_config_file"
        sed -i "s|databases_root = \"/user/home/.metagear/databases\"|databases_root = \"$HOME/.metagear/databases\"|g" "$user_config_file"

        cp $1/templates/metagear.env $user_env_file

        echo ""
        echo "It seems this is the first timeMetaGEAR.."
        echo ""
        check_requirements
        echo ""
        echo "   - User configuration was created in ~/.metagear/metagear.config"
        echo "   - Environment file was created in ~/.metagear/metagear.env"
        echo ""
        echo "IMPORTANT: Please review these file before re-launching the MetaGEAR pipeline."
        echo ""

        check_requirements

        exit 0

    fi

}


detect_container_tool() {
    # Check for singularity first and return it if found.
    if command -v singularity >/dev/null 2>&1; then
        echo "singularity"
    # Else check for docker and return it if found.
    elif command -v docker >/dev/null 2>&1; then
        echo "docker"
    # Return "none" if neither is available.
    else
        echo "Error: MetaGEAR requires Singularity (recommended) or Docker. Please install one of those in your system." >&2
        echo "  Singularity: https://docs.sylabs.io/guides/3.0/user-guide/installation.html" >&2
        echo "  Docker: https://docs.docker.com/engine/install/" >&2
        exit 1
    fi
}

