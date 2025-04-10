#!/bin/bash
# merge_configs.sh
# Usage: ./merge_configs.sh file1 file2 [file3 ...]
# Merges configuration files in order; later file values/blocks override earlier ones.
# This version supports three sections: params, profiles, and process.

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 file1 file2 [file3 ...]"
    exit 1
fi

# Declare associative arrays for each section
declare -A params_map    # For key = value assignments inside the params section
declare -A proc_map      # For withName blocks inside the process section
declare -A profiles_map  # For profile blocks inside the profiles section

# Process files in the order provided.
for file in "$@"; do
    if [ ! -f "$file" ]; then
        echo "Error: file '$file' not found." >&2
        exit 1
    fi

    # current_section will be either "params", "process", or "profiles"
    current_section=""
    while IFS= read -r line; do
        # Remove leading and trailing whitespace
        trimmed=$(echo "$line" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')

        # Detect section headers
        if [[ "$trimmed" == "params {" ]]; then
            current_section="params"
            continue
        elif [[ "$trimmed" == "process {" ]]; then
            current_section="process"
            continue
        elif [[ "$trimmed" == "profiles {" ]]; then
            current_section="profiles"
            continue
        fi

        # Detect a top-level closing "}"
        if [[ "$trimmed" == "}" ]]; then
            if [[ "$current_section" == "params" || "$current_section" == "process" || "$current_section" == "profiles" ]]; then
                current_section=""
            fi
            continue
        fi

        # Handle each section accordingly
        if [[ "$current_section" == "params" ]]; then
            # Skip blank lines or comments
            if [[ -z "$trimmed" || "$trimmed" =~ ^/\* ]]; then
                continue
            fi
            # Match assignments of the form: key = value
            if [[ $trimmed =~ ^([a-zA-Z_][a-zA-Z0-9_]*)[[:space:]]*=[[:space:]]*(.*)$ ]]; then
                key="${BASH_REMATCH[1]}"
                value="${BASH_REMATCH[2]}"
                params_map["$key"]="$value"
            fi

        elif [[ "$current_section" == "process" ]]; then
            # Look for process sub-blocks that start with "withName:"
            if [[ "$trimmed" =~ ^withName: ]]; then
                if [[ "$trimmed" =~ ^withName:[[:space:]]*['"]?([^'":]+)['"]?[[:space:]]*\{ ]]; then
                    block_name="${BASH_REMATCH[1]}"
                else
                    continue
                fi
                block_content="$line"$'\n'
                # Count braces in the current line.
                open_braces=$(grep -o "{" <<< "$line" | wc -l)
                close_braces=$(grep -o "}" <<< "$line" | wc -l)
                brace_count=$((open_braces - close_braces))
                # Continue reading until the block is balanced.
                while [ $brace_count -gt 0 ] && IFS= read -r next_line; do
                    block_content+="$next_line"$'\n'
                    open_braces=$(grep -o "{" <<< "$next_line" | wc -l)
                    close_braces=$(grep -o "}" <<< "$next_line" | wc -l)
                    brace_count=$((brace_count + open_braces - close_braces))
                done
                # Override any previous block with the same name.
                proc_map["$block_name"]="$block_content"
            fi

        elif [[ "$current_section" == "profiles" ]]; then
            # In profiles section, each profile block starts with a line like:
            #   profile_name {
            if [[ "$trimmed" =~ ^([a-zA-Z0-9_]+)[[:space:]]*\{ ]]; then
                profile_name="${BASH_REMATCH[1]}"
                block_content="$line"$'\n'
                open_braces=$(grep -o "{" <<< "$line" | wc -l)
                close_braces=$(grep -o "}" <<< "$line" | wc -l)
                brace_count=$((open_braces - close_braces))
                while [ $brace_count -gt 0 ] && IFS= read -r next_line; do
                    block_content+="$next_line"$'\n'
                    open_braces=$(grep -o "{" <<< "$next_line" | wc -l)
                    close_braces=$(grep -o "}" <<< "$next_line" | wc -l)
                    brace_count=$((brace_count + open_braces - close_braces))
                done
                profiles_map["$profile_name"]="$block_content"
            fi
        fi
    done < "$file"
done

# Write the merged configuration to merged.txt
output_file="metagear_run.config"
{
    echo "params {"
    for key in "${!params_map[@]}"; do
        echo "    $key = ${params_map[$key]}"
    done
    echo "}"
    echo ""
    echo "profiles {"
    for key in "${!profiles_map[@]}"; do
        echo "${profiles_map[$key]}"
    done
    echo "}"
    echo ""
    echo "process {"
    for key in "${!proc_map[@]}"; do
        echo "${proc_map[$key]}"
    done
    echo "}"
} > "$output_file"

echo "Merged configuration saved to $output_file"
