#!/bin/bash

# Define the source directory for scripts
scripts_dir="./scripts"

# Check if the scripts directory exists
if [ ! -d "$scripts_dir" ]; then
    echo "Error: Scripts directory '$scripts_dir' not found."
    exit 1
fi

# Define the list of files to copy
files_to_copy=("mat_generation.jl" "matGain_generation.jl" "rl.jl" "rl_gold.jl" "rl.sh" "rl_gold.sh")

# Verify that all files exist in the scripts directory
for file in "${files_to_copy[@]}"; do
    if [ ! -f "$scripts_dir/$file" ]; then
        echo "Error: File '$file' not found in '$scripts_dir'."
        exit 1
    fi
done

# Check if input_data.txt exists
input_data_file="$scripts_dir/input_data.txt"
if [ ! -f "$input_data_file" ]; then
    echo "Error: File 'input_data.txt' not found in '$scripts_dir'."
    exit 1
fi

# Ensure all .sh scripts are executable
for file in "${files_to_copy[@]}"; do
    if [[ "$file" == *.sh ]]; then
        chmod +x "$scripts_dir/$file"
    fi
done

# Iterate over all sample_* directories
for dir in sample_*; do
    if [ -d "$dir" ]; then
        echo "Copying files to directory: $dir"
        # Copy each script to the target directory
        for file in "${files_to_copy[@]}"; do
            cp "$scripts_dir/$file" "$dir"
            # Ensure .sh files are executable in the target directory
            if [[ "$file" == *.sh ]]; then
                chmod +x "$dir/$file"
            fi
        done
        # Copy input_data.txt to the InputData subdirectory
        input_data_dir="$dir/InputData"
        mkdir -p "$input_data_dir" # Create InputData directory if it doesn't exist
        cp "$input_data_file" "$input_data_dir"
    fi
done

echo "All files have been copied and updated in the sample directories."

