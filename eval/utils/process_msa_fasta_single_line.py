import argparse
import os

def process_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        sequence = ''
        for line in infile:
            if line.startswith('>'):
                if sequence:
                    outfile.write(sequence + '\n')
                outfile.write(line)
                sequence = ''
            else:
                # capitalize sequence
                sequence += line.strip().upper()
        if sequence:
            outfile.write(sequence + '\n')

def main():
    parser = argparse.ArgumentParser(description='Convert multi-line FASTA to single-line FASTA')
    parser.add_argument('-i', '--input-file', type=str, help='Input multi-line FASTA file')
    parser.add_argument('-o', '--output-file', type=str, help='Output single-line FASTA file')
    parser.add_argument('-d', '--dir', type=str, help='Directory containing multi-line FASTA files')
    
    args = parser.parse_args()

    if args.dir:
        for file in os.listdir(args.dir):
            if file.endswith('.fasta'):
                input_file = os.path.join(args.dir, file)
                out_dir = os.path.join(args.dir, 'processed')
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)
                output_file = os.path.join(out_dir, file)
                process_fasta(input_file, output_file)
                print(f'Processed {input_file}')
    else:
        process_fasta(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
