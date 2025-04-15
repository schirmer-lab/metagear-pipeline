#!/bin/bash
# merge_configuration.sh (Bash 3 Compatible)
# Usage: ./merge_configuration.sh file1 file2 [file3 ...]
# Merges configuration files in order; later file values/blocks override earlier ones.
# Supports three sections: params, profiles, and process.

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 file1 file2 [file3 ...]" >&2
    exit 1
fi

# Initialize parallel arrays for each section
params_keys=()
params_values=()

proc_keys=()
proc_values=()

profiles_keys=()
profiles_values=()

# Function to set or update a key/value pair for params
set_param() {
    local key="$1"
    local value="$2"
    local i
    for i in "${!params_keys[@]}"; do
        if [ "${params_keys[$i]}" = "$key" ]; then
            params_values[$i]="$value"
            return
        fi
    done
    params_keys+=("$key")
    params_values+=("$value")
}

# Function to set or update a process block
set_proc() {
    local key="$1"
    local value="$2"
    local i
    for i in "${!proc_keys[@]}"; do
        if [ "${proc_keys[$i]}" = "$key" ]; then
            proc_values[$i]="$value"
            return
        fi
    done
    proc_keys+=("$key")
    proc_values+=("$value")
}

# Function to set or update a profile block
set_profile() {
    local key="$1"
    local value="$2"
    local i
    for i in "${!profiles_keys[@]}"; do
        if [ "${profiles_keys[$i]}" = "$key" ]; then
            profiles_values[$i]="$value"
            return
        fi
    done
    profiles_keys+=("$key")
    profiles_values+=("$value")
}

# Process input files in order
for file in "$@"; do
    if [ ! -f "$file" ]; then
        echo "Warning: file '$file' not found. Skipping..." >&2
        continue
    fi

    current_section=""
    while IFS= read -r line; do
        # Trim leading and trailing whitespace
        trimmed=$(echo "$line" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')

        # Detect section headers
        case "$trimmed" in
            "params {")
                current_section="params"
                continue
                ;;
            "process {")
                current_section="process"
                continue
                ;;
            "profiles {")
                current_section="profiles"
                continue
                ;;
            "}")
                if [ "$current_section" = "params" ] || [ "$current_section" = "process" ] || [ "$current_section" = "profiles" ]; then
                    current_section=""
                fi
                continue
                ;;
        esac

        # Process each section accordingly
        if [ "$current_section" = "params" ]; then
            # Skip blank lines or comments (assumed to start with '/')
            case "$trimmed" in
                "") continue ;;
                /*) continue ;;
            esac
            # Match assignments of the form: key = value
            if echo "$trimmed" | grep -qE '^[a-zA-Z_][a-zA-Z0-9_]*[[:space:]]*=[[:space:]]*'; then
                key=$(echo "$trimmed" | sed -E 's/^([a-zA-Z_][a-zA-Z0-9_]*)[[:space:]]*=.*/\1/')
                value=$(echo "$trimmed" | sed -E 's/^[a-zA-Z_][a-zA-Z0-9_]*[[:space:]]*=[[:space:]]*(.*)/\1/')
                set_param "$key" "$value"
            fi

        elif [ "$current_section" = "process" ]; then
            if echo "$trimmed" | grep -q '^withName:'; then
                # Extract block name using sed; supports optional quotes
                block_name=$(echo "$trimmed" | sed -E "s/^withName:[[:space:]]*['\"]?([^'\"]+)['\"]?[[:space:]]*\{.*/\1/")
                block_content="$line"$'\n'
                open_braces=$(echo "$line" | grep -o "{" | wc -l)
                close_braces=$(echo "$line" | grep -o "}" | wc -l)
                brace_count=$((open_braces - close_braces))
                # Accumulate lines until braces balance
                while [ $brace_count -gt 0 ] && IFS= read -r next_line; do
                    block_content="$block_content$next_line"$'\n'
                    open_braces=$(echo "$next_line" | grep -o "{" | wc -l)
                    close_braces=$(echo "$next_line" | grep -o "}" | wc -l)
                    brace_count=$((brace_count + open_braces - close_braces))
                done
                set_proc "$block_name" "$block_content"
            fi

        elif [ "$current_section" = "profiles" ]; then
            if echo "$trimmed" | grep -qE '^[a-zA-Z0-9_]+[[:space:]]*\{'; then
                profile_name=$(echo "$trimmed" | sed -E 's/^([a-zA-Z0-9_]+)[[:space:]]*\{.*/\1/')
                block_content="$line"$'\n'
                open_braces=$(echo "$line" | grep -o "{" | wc -l)
                close_braces=$(echo "$line" | grep -o "}" | wc -l)
                brace_count=$((open_braces - close_braces))
                while [ $brace_count -gt 0 ] && IFS= read -r next_line; do
                    block_content="$block_content$next_line"$'\n'
                    open_braces=$(echo "$next_line" | grep -o "{" | wc -l)
                    close_braces=$(echo "$next_line" | grep -o "}" | wc -l)
                    brace_count=$((brace_count + open_braces - close_braces))
                done
                set_profile "$profile_name" "$block_content"
            fi
        fi
    done < "$file"
done

# Write the merged configuration to the output file

{
    echo "params {"
    for i in "${!params_keys[@]}"; do
        echo "    ${params_keys[$i]} = ${params_values[$i]}"
    done
    echo "}"
    echo ""
    echo "profiles {"
    for i in "${!profiles_keys[@]}"; do
        echo "${profiles_values[$i]}"
    done
    echo "}"
    echo ""
    echo "process {"
    for i in "${!proc_keys[@]}"; do
        echo "${proc_values[$i]}"
    done
    echo "}"
}

