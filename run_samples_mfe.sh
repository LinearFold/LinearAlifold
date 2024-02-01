#!/bin/bash

# Use the first command line argument as the input directory
input_directory="$1"
output_directory="$2"
em_value="$3"  # the -e argument value (1 or 2)
multi_approx="$4"
cutoff="$5"
beta="$6"
delta="$7"
partition_mode="$8"

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

        # check if partition mode is enabled
        if [ "$partition_mode" = 1 ] ; then
            # Run LinearAlifold with partition mode enabled
            cat "$file" | python3 ./linearalifold.py --verbose --em $em_value --multi-approx $multi_approx -ct $cutoff -bt $beta -dt $delta --partition > "$output_file"
        else
            # Run LinearAlifold with partition mode disabled
            cat "$file" | python3 ./linearalifold.py --verbose --em $em_value --multi-approx $multi_approx -ct $cutoff -bt $beta -dt $delta > "$output_file"
        fi

        # Record the end time
        end_time=$(date +%s)

        # Calculate the runtime
        runtime=$((end_time - start_time))
        echo "Processed $file in $runtime seconds."
    fi
done

echo "Processing complete."
