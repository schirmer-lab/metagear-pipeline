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
        echo "-> It seems this is the first timeMetaGEAR is running in this system."
        echo "   - User configuration was created in ~/.metagear/metagear.config"
        echo "   - Environment file was created in ~/.metagear/metagear.env"
        echo ""
        echo "Please review these file before re-launching the MetaGEAR pipeline."
        echo ""
        exit 0

    fi

}


