#!/bin/bash

# This script should be called with two arguments.
# The first argument is the directory containing the input files.
# The second argument is the value for the -e (energy model) option [should be 1 (Original) or 2 (LinearAliFold)].
# Example usage: ./run_samples.sh my_input_directory 1

# Use the first command line argument as the input directory
input_directory="$1"
output_directory="$2"

# Create a directory for the output files
mkdir -p "$output_directory"

echo "Results will be stored in $output_directory"

for file in "$input_directory"/*
do
    if [ -f "$file" ]; then
        # Get the filename without the directory
        filename=$(basename -- "$file")
        # Create a new output file for each input file in the output directory
        output_file="${output_directory}/${filename}/"
        echo "Processing $file..."
        
        # Record the start time
        start_time=$(date +%s)

        ./linearturbofold -i "$file" -o "$output_file"

        # Record the end time
        end_time=$(date +%s)

        # Calculate the runtime
        runtime=$((end_time - start_time))
        echo "Processed $file in $runtime seconds."
    fi
done

echo "Processing complete."
