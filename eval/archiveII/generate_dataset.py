import os
from pathlib import Path

# Define a function to parse the FASTA file and write the sequence and structure into separate files
def parse_and_write_fasta(file_path, sequences_dir, structures_dir):
    with open(file_path, 'r') as file:
        while True:
            # Read the sequence name line
            name_line = file.readline().strip()
            if not name_line:  # If line is empty, we are done
                break
            
            # Extract the sequence name
            seq_name = name_line.split('/')[1]
            
            # Read the sequence line
            sequence = file.readline().strip()
            
            # Read the structure line
            structure = file.readline().strip()

            # Write the sequence to its file
            sequence_file = sequences_dir / f"{seq_name}.txt"
            with open(sequence_file, 'w') as sf:
                sf.write(">" + seq_name + "\n" + sequence + "\n")

            # Write the structure to its file, if it exists
            if structure:
                structure_file = structures_dir / f"{seq_name}.txt"
                with open(structure_file, 'w') as stf:
                    stf.write(">" + seq_name + "\n" + structure + "\n")

# Path to the input FASTA file
file_path = './data_raw.fasta'

# Define the sequences and structures directories
sequences_dir = Path('./sequences')
structures_dir = Path('./structures')

# Create the directories if they don't already exist
sequences_dir.mkdir(parents=True, exist_ok=True)
structures_dir.mkdir(parents=True, exist_ok=True)

# Call the function to parse the file and create sequence and structure files
parse_and_write_fasta(file_path, sequences_dir, structures_dir)
