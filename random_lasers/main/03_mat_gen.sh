#!/bin/bash

# Check if Julia is installed
if ! command -v julia &> /dev/null; then
    echo "Julia is not installed or not in the PATH."
    exit 1
fi

# Find all directories matching the pattern sample_*
find . -type d -name "sample_*" | parallel --bar '
    echo "Processing {}"
    cd "{}" || exit
    julia mat_generation.jl
    julia matGain_generation.jl
'
