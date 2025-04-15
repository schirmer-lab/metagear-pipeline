#!/bin/bash

# Support vairables
script_dir="$(cd "$(dirname "$0")" && pwd)"
launch_dir="$PWD"

${script_dir}/.commander/main.sh "$@"
