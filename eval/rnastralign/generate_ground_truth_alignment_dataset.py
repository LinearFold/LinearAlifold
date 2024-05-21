import os
import argparse

from find_ground_truth_alignment import search_fasta_files

def copy_and_query_fasta(input_dir, output_dir, database_dir):
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Walk through the input directory to find all .fasta files
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.fasta'):
                # if file name starts with 23S skip
                if file.startswith("23S"):
                    continue

                # Determine the full input file path
                input_file_path = os.path.join(root, file)

                # Determine the corresponding output file path
                relative_path = os.path.relpath(root, input_dir)
                output_file_dir = os.path.join(output_dir, relative_path)
                if not os.path.exists(output_file_dir):
                    os.makedirs(output_file_dir)

                output_file_path = os.path.join(output_file_dir, file)

                # Read from input file and write to output file
                with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
                    lines = input_file.readlines()
                    print("Processing", input_file_path)
                    for i in range(len(lines)):
                        # if line is empty or doesn't start with '>' skip
                        if not lines[i].strip() or not lines[i].startswith(">"):
                            continue
                        query = lines[i].strip().split('/')[-1]
                        query = ".".join(query.split('.')[:-1])
                        result = search_fasta_files(database_dir, query)

                        # print(f"Query: {query}, Matched Line: {result[0][0]}")

                        output_file.write(lines[i].strip() + "\n")
                        output_file.write(result[0][1] + "\n")

def main():
    parser = argparse.ArgumentParser(description="Extract query from each line in .fasta files and write the next line if it contains the query.")
    parser.add_argument('--database', help="Directory to search for .fasta files", default='./../database')
    parser.add_argument('--input', help="Input directory to search for .fasta files", default='./aln')
    parser.add_argument('--output', help="Output directory where the resulting files will be created")
    
    args = parser.parse_args()
    
    copy_and_query_fasta(args.input, args.output, args.database)

if __name__ == "__main__":
    main()
