#!/bin/bash

# Use the first command line argument as the input directory
input_directory="$1"
output_directory="$2"
em_value="$3"  # the -e argument value (1 or 2)
cutoff="$4"
beta="$5"
delta="$6"

# Create a directory for the output files
mkdir -p "$output_directory"

echo "Results will be stored in $output_directory"

for file in "$input_directory"/*
do
    if [ -f "$file" ]; then
        # Get the filename without the directory
        filename=$(basename -- "$file")
        # Create a new output file for each input file in the output directory
        output_file="${output_directory}/${filename}.txt"
        echo "Processing $file..."

        # Record the start time
        start_time=$(date +%s)

        cat "$file" | python3 ./linearalifold.py --verbose -e $em_value -c $cutoff -y $beta -z $delta > "$output_file"

        # Record the end time
        end_time=$(date +%s)

        # Calculate the runtime
        runtime=$((end_time - start_time))
        echo "Processed $file in $runtime seconds."
    fi
done

echo "Processing complete."
