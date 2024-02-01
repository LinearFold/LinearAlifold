#!/bin/bash

# This script should be called with two arguments.
# The first argument is the directory containing the input files.
# The second argument is the value for the -e (energy model) option [should be 1 (Original) or 2 (LinearAliFold)].
# Example usage: ./run_samples.sh my_input_directory 1

# Use the first command line argument as the input directory
input_directory="$1"
output_directory="$2"
em_value="$3"  # the -e argument value (1 or 2)
multi_approx="$4"
cutoff="$5"
beta="$6"
delta="$7"

# Create a directory for the output files
mkdir -p "$output_directory"

# create log, mea, bpp folders in output directory
mkdir -p "$output_directory/log"
mkdir -p "$output_directory/mea"
mkdir -p "$output_directory/bpp"
mkdir -p "$output_directory/threshknot"

echo "Results will be stored in $output_directory"

for file in "$input_directory"/*
do
    if [ -f "$file" ]; then
        # Get the filename without the directory
        filename=$(basename -- "$file")
        # Create a new output file for each input file in the output directory
        log_file="${output_directory}/log/${filename}.txt"
        mea_file="${output_directory}/mea/${filename}.txt"
        bpp_file="${output_directory}/bpp/${filename}.txt"
        threshknot_file="${output_directory}/threshknot/${filename}.txt"

        echo "Processing $file..."
        
        # Record the start time
        start_time=$(date +%s)

        cat "$file" | ./linearalifold.py --partition --em $em_value --bpp-file $bpp_file --mea-file $mea_file --threshknot-file $threshknot_file\
         --multi-approx $multi_approx -ct $cutoff -bt $beta -dt $delta --verbose > $log_file

        # Record the end time
        end_time=$(date +%s)

        # Calculate the runtime
        runtime=$((end_time - start_time))
        echo "Processed $file in $runtime seconds."
    fi
done

echo "Processing complete."
