import os
import argparse

def search_fasta_files(directory, query):
    results = []
    # Walk through all directories and files recursively
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.fasta'):
                # Construct the full file path
                filepath = os.path.join(root, file)
                with open(filepath, 'r') as fasta_file:
                    lines = fasta_file.readlines()
                    # Iterate through the lines in the file
                    for i in range(len(lines) - 1):
                        if query in lines[i]:
                            # Append the matched line and the next line
                            results.append((lines[i].strip(), lines[i+1].strip()))
                            return results
    
    raise ValueError(f"Query '{query}' not found in any .fasta files in directory '{directory}'")

def main():
    parser = argparse.ArgumentParser(description="Search for lines containing a query in FASTA files and extract the next line.")
    parser.add_argument('--directory', help="Directory to search for .fasta files", default='./../database')
    parser.add_argument('--query', help="Query string to search for within the .fasta files")
    
    args = parser.parse_args()
    
    results = search_fasta_files(args.directory, args.query)
    for result in results:
        print(f"Matched Line: {result[0]}\nNext Line: {result[1]}\n")

if __name__ == "__main__":
    main()
