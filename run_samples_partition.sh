#!/bin/bash

# This script should be called with two arguments.
# The first argument is the directory containing the input files.
# The second argument is the value for the -e (energy model) option [should be 1 (Original) or 2 (LinearAliFold)].
# Example usage: ./run_samples.sh my_input_directory 1

# Use the first command line argument as the input directory
input_directory="$1"
output_directory="$2"
em_value="$3"  # the -e argument value (1 or 2)
cutoff="$4"
beta="$5"
delta="$6"

# Create a directory for the output files
mkdir -p "$output_directory"

# create log, mea, bpp folders in output directory
mkdir -p "$output_directory/log"
mkdir -p "$output_directory/mea"
mkdir -p "$output_directory/bpp"
mkdir -p "$output_directory/threshknot"
mkdir -p "$output_directory/centroid"

echo "Results will be stored in $output_directory"

for file in $(echo "$input_directory"/* | tr ' ' '\n' | sort -V);
do
    if [ -f "$file" ]; then
        # Get the filename without the directory
        filename=$(basename -- "$file")
        # Create a new output file for each input file in the output directory
        log_file="${output_directory}/log/${filename}.txt"
        mea_file="${output_directory}/mea/${filename}.txt"
        bpp_file="${output_directory}/bpp/${filename}.txt"
        threshknot_file="${output_directory}/threshknot/${filename}.txt"
        centroid_file="${output_directory}/centroid/${filename}.txt"

        echo "Processing $file..."
        
        # Record the start time
        start_time=$(date +%s)

        cat "$file" | ./linearalifold.py --partition -e $em_value -B $bpp_file -M $mea_file -T $threshknot_file\
         -C $centroid_file -c $cutoff -y $beta -z $delta --verbose > $log_file

        # Record the end time
        end_time=$(date +%s)

        # Calculate the runtime
        runtime=$((end_time - start_time))
        echo "Processed $file in $runtime seconds."
    fi
done

echo "Processing complete."
