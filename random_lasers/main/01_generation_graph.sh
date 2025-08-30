#!/bin/bash

# Prompt the user for the range of samples to generate
read -p "Enter the range of samples to generate (e.g., 9-200): " range_input

# Validate the input format
if ! [[ "$range_input" =~ ^[0-9]+-[0-9]+$ ]]; then
    echo "Error: Please enter a valid range in the format start-end (e.g., 9-200)."
    exit 1
fi

# Extract the start and end of the range
start_sample=$(echo "$range_input" | cut -d'-' -f1)
end_sample=$(echo "$range_input" | cut -d'-' -f2)

# Validate that start and end are positive integers and that start <= end
if ! [[ "$start_sample" =~ ^[0-9]+$ ]] || ! [[ "$end_sample" =~ ^[0-9]+$ ]] || [ "$start_sample" -le 0 ] || [ "$end_sample" -le 0 ] || [ "$start_sample" -gt "$end_sample" ]; then
    echo "Error: Invalid range. Ensure start and end are positive integers, and start <= end."
    exit 1
fi

# Define the base directory and Julia script path
base_dir=$(pwd)
julia_script="$base_dir/scripts/generation_graph.jl"  # Script Julia nella directory scripts

# Debug: Print paths
echo "Using Julia script: $julia_script"

# Check if the Julia script exists
if [ ! -f "$julia_script" ]; then
    echo "Error: Julia script not found at $julia_script."
    exit 1
fi

# Generate the list of sample directories
sample_dirs=$(seq -f "sample_%03g" $start_sample $end_sample)

# Debug: Print new sample directories to be created
echo "Sample directories to create: $sample_dirs"

# Create new sample directories
for dir in $sample_dirs; do
    mkdir -p "$base_dir/$dir"
done

# Record the start time
start_time=$(date +%s)

# Execute the Julia script for each new directory sequentially
for dir in $sample_dirs; do
    sample_path="$base_dir/$dir"
    echo "Running Julia script for directory: $sample_path"
    julia "$julia_script" "$sample_path"
done

# Record the end time
end_time=$(date +%s)

# Calculate and display the total execution time
total_time=$((end_time - start_time))
echo "All samples from $start_sample to $end_sample have been generated!"
echo "Total execution time: $total_time seconds"

