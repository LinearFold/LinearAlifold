#!/bin/bash

input_directory="$1"
output_directory="$2"
run_part="$3"       # 1 to run centroid, 0 to run mfe

# Create a directory for the output files
mkdir -p "$output_directory"
mkdir -p "$output_directory/output"
mkdir -p "$output_directory/log"

echo "Results will be stored in $output_directory"

# Create an array of sorted file names using version sort
readarray -t sorted_files < <(find "$input_directory" -maxdepth 1 -type f | sort -V)

for file in "${sorted_files[@]}"
do
    if [ -f "$file" ]; then
        # Get the filename without the directory
        filename=$(basename -- "$file").txt

        echo "Processing $file..."
        
        # Record the start time
        start_time=$(date +%s)

        if [ "$run_part" = 1 ] ; then
            python3 run.py -i "$file" -o "$output_directory/output/$filename" --pt > "$output_directory/log/$filename"
        else
            python3 run.py -i "$file" -o "$output_directory/output/$filename" > "$output_directory/log/$filename"
        fi

        # Record the end time
        end_time=$(date +%s)

        # Calculate the runtime
        runtime=$((end_time - start_time))
        echo "Processed $file in $runtime seconds."
    fi
done

echo "Processing complete."
