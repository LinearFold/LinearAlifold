import os
import argparse

def remove_gaps_in_sequences(directory):
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path):
            # Read the original content of the file
            with open(file_path, 'r') as file:
                lines = file.readlines()

            # Process each line and remove gaps if the line is not a FASTA header
            with open(file_path, 'w') as file:
                for line in lines:
                    if not line.startswith(">"):
                        line = line.replace("-", "")  # Remove gaps
                    file.write(line)
            print(f"Updated file: {file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove gaps from sequence lines in all files of a directory.")
    parser.add_argument("directory", help="Path to the directory containing the files")

    args = parser.parse_args()

    remove_gaps_in_sequences(args.directory)
