#!/usr/bin/env bash

# Function to get total installed memory in GB.
# It reads /proc/meminfo for MemTotal (in kB) and converts it to GB.
get_total_memory_gb() {
    local mem_total_kb
    mem_total_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    # Divide by 1048576 (1024*1024) to convert kB to GB.
    local mem_total_gb
    mem_total_gb=$(awk -v mem="$mem_total_kb" 'BEGIN {printf "%.2f", mem/1048576}')
    echo "$mem_total_gb"
}

# Function to get the number of available CPUs.
# It uses the nproc command, which is available by default on Ubuntu.
get_cpu_count() {
    if command -v nproc >/dev/null 2>&1; then
        nproc
    else
        # Fallback: if nproc is not available, default to 1.
        echo "1"
    fi
}
