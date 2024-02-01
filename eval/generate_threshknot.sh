input_directory="$1"
output_directory="$2"
threshold="$3"

# if output directory does not exist, create it
if [ ! -d "$output_directory" ]; then
    mkdir -p "$output_directory"
fi

for file in "$input_directory"/*
do
    if [ -f "$file" ]; then
        # Get the filename without the directory
        filename=$(basename -- "$file")
        # Create a new output file for each input file in the output directory
        threshknot_file="${output_directory}/${filename}"

        echo "Processing $file..."

        python3 generate_threshknot.py -i $file -o $threshknot_file --threshold $threshold
    fi
done