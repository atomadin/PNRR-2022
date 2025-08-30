#!/bin/bash

# Specify the numbers directly in the script, example: target_input="2-5,7"
target_input="1-100"

# Function to parse ranges and comma-separated values
parse_input() {
    local input="$1"
    local -n result=$2
    result=()  # Initialize the array

    # Split the input into parts separated by commas
    IFS=',' read -ra parts <<< "$input"
    for part in "${parts[@]}"; do
        if [[ "$part" =~ ^[0-9]+-[0-9]+$ ]]; then
            # Expand ranges (e.g., 2-5) into individual numbers
            IFS='-' read -r start end <<< "$part"
            for ((i = start; i <= end; i++)); do
                result+=("$i")
            done
        elif [[ "$part" =~ ^[0-9]+$ ]]; then
            # Add single numbers directly
            result+=("$part")
        else
            # Handle invalid input
            echo "Invalid input: $part"
            exit 1
        fi
    done
}

# Parse the input to generate the target_n array
target_n=()
parse_input "$target_input" target_n

# Print the parsed numbers
echo "Target n values: ${target_n[*]}"

# Iterate over the parsed numbers
for n in "${target_n[@]}"; do
    folder=$(printf "sample_%03d" "$n")  # Format folder name with leading zeros
    script="rl_gold.sh"                      # Script name to execute

    if [[ -d "$folder" ]]; then
        # Enter the folder if it exists
        echo "Entering folder $folder..."
        cd "$folder" || exit 1

        if [[ -x "$script" ]]; then
            # Run the script if it is executable
            echo "Running $script in folder $folder..."
            ./"$script"
        else
            # Notify if the script is not executable or missing
            echo "Script $script is not executable or does not exist in folder $folder."
        fi

        # Return to the main directory
        cd ..
    else
        # Notify if the folder does not exist
        echo "Folder $folder does not exist."
    fi
done

# Confirm the operation is completed
echo "Operation completed."
