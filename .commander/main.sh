#!/usr/bin/env bash
# .commander/main.sh

# Resolve script directory and source common functions
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
LAUNCH_DIR="$PWD"
source "${SCRIPT_DIR}/lib/common.sh"
source "${SCRIPT_DIR}/lib/workflows.sh"

# Ensure a command is provided
if [ $# -eq 0 ]; then
    usage
fi

COMMAND="$1"
shift

check_command "$COMMAND"
check_metagear_home $SCRIPT_DIR

mkdir -p $LAUNCH_DIR/.metagear

custom_config_files=( $SCRIPT_DIR/../conf/metagear/$COMMAND.config $HOME/.metagear/metagear.config )
metagear_config_files=( $SCRIPT_DIR/../conf/metagear/*.config )
all_config_files=( "${metagear_config_files[@]}" "${custom_config_files[@]}" )

$SCRIPT_DIR/lib/merge_configuration.sh ${all_config_files[@]} > $LAUNCH_DIR/.metagear/$COMMAND.config

nf_cmd_workflow_part=$(run_workflows $COMMAND $@)

cat $HOME/.metagear/metagear.env > $LAUNCH_DIR/metagear_$COMMAND.sh

echo "" >> $LAUNCH_DIR/metagear_$COMMAND.sh
echo "nextflow run $SCRIPT_DIR/../main.nf \\
        $nf_cmd_workflow_part \\
        -c $LAUNCH_DIR/.metagear/$COMMAND.config \\
        \$RUN_PROFILES -w \\
        \$NF_WORK -resume" >> $LAUNCH_DIR/metagear_$COMMAND.sh
echo "" >> $LAUNCH_DIR/metagear_$COMMAND.sh

chmod +x $LAUNCH_DIR/metagear_$COMMAND.sh

$LAUNCH_DIR/metagear_$COMMAND.sh
