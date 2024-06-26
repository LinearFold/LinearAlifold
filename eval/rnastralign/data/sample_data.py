import os
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import json

def is_rna_sequence(line):
    """Check if a line contains a valid RNA sequence, ignoring the trailing '1'."""
    valid_nucleotides = set("AUCG")
    cleaned_line = line.strip().rstrip('1')
    return all(char in valid_nucleotides for char in cleaned_line)

def sample_rna_sequences(directory, k, family_name, min_length=1):
    """
    Sample a specified number of RNA sequences from a directory containing .seq files
    that are greater than a given minimum length.
    
    Args:
    directory (str): Path to the directory containing .seq files.
    k (int): Number of sequences to sample for each MSA.
    min_length (int): Minimum length of sequences to sample.
    
    Returns:
    list: List of sampled SeqRecord objects.
    """
    seq_files = [f for f in os.listdir(directory) if f.endswith('.seq')]

    if family_name == 'SRP':
        min_length = 200
    if family_name == '16S':
        min_length = 1200
    
    valid_sequences = []
    
    for file in seq_files:
        with open(os.path.join(directory, file), 'r') as f:
            lines = f.readlines()
            sequence = None
            for line in lines:
                if is_rna_sequence(line):
                    sequence = line.strip().rstrip('1')
                    break
            if sequence and len(sequence) >= min_length:
                seq_id = os.path.join(directory, file)
                if seq_id.split('.')[-1].strip() == 'seq':
                    seq_id = seq_id[:-4]
                record = SeqRecord(
                    Seq(sequence),
                    id=seq_id + '.ct',
                    description=""
                )
                valid_sequences.append(record)
    
    if len(valid_sequences) < k:
        raise ValueError(f"Not enough sequences longer than {min_length} found in the directory")
    
    sampled_sequences = random.sample(valid_sequences, k)
    
    return sampled_sequences

def create_msa_fasta(sequences, output_file):
    """
    Create a MSA in FASTA format from a list of sequences.
    
    Args:
    sequences (list): List of SeqRecord objects.
    output_file (str): Path to the output FASTA file.
    """
    with open(output_file, 'w') as f:
        for seq_record in sequences:
            f.write(f">{seq_record.id}\n{str(seq_record.seq)}\n")

def align_with_mafft(input_file, output_file):
    """
    Align sequences in the input FASTA file using MAFFT and save to the output file.
    
    Args:
    input_file (str): Path to the input FASTA file.
    output_file (str): Path to the output aligned FASTA file.
    """
    result = subprocess.run(['mafft', '--auto', input_file], stdout=open(output_file, 'w'))
    # subprocess.run waits for the command to complete by default
    if result.returncode != 0:
        raise RuntimeError(f"MAFFT alignment failed for {input_file}. Return code: {result.returncode}")


def format_fasta_single_line(input_file):
    """
    Format sequences in a FASTA file to be in a single line and capitalized.
    
    Args:
    input_file (str): Path to the input FASTA file.
    """
    records = []
    for record in SeqIO.parse(input_file, "fasta"):
        record.seq = Seq(str(record.seq).upper())
        records.append(record)
    
    with open(input_file, 'w') as f:
        for record in records:
            f.write(f">{record.id}\n{str(record.seq)}\n")

def generate_msas(input_dir, output_dir, k, count, dataset_name, family_name, subfamily_name, align):
    """
    Generate multiple MSAs in FASTA format from sampled RNA sequences.
    
    Args:
    input_dir (str): Path to the directory containing .seq files.
    output_dir (str): Path to the directory where output files will be saved.
    k (int): Number of sequences per MSA.
    count (int): Number of MSAs to generate.
    dataset_name (str): Name of the dataset.
    family_name (str): Name of the family.
    subfamily_name (str): Name of the subfamily to include as a prefix in sequence IDs.
    align (int): Alignment option (0: no alignment, 1: do alignment, 2: generate both versions).
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if align == 2:
        aln_dir = os.path.join(output_dir, "aln")
        no_aln_dir = os.path.join(output_dir, "no_aln")
        os.makedirs(aln_dir, exist_ok=True)
        os.makedirs(no_aln_dir, exist_ok=True)
    else:
        aln_dir = output_dir
        no_aln_dir = output_dir
    
    for i in range(1, count + 1):
        base_filename = f"{family_name}.{dataset_name}.k_{k}_{i}.fasta"
        
        if align == 0:
            output_file = os.path.join(no_aln_dir, base_filename)
            sampled_sequences = sample_rna_sequences(input_dir, k, family_name)
            create_msa_fasta(sampled_sequences, output_file)
            print(f"Non-aligned MSA file created: {output_file}")

        elif align == 1:
            aligned_filename = f"{family_name}.{dataset_name}.k_{k}_{i}.aln.fasta"
            aligned_output_file = os.path.join(aln_dir, aligned_filename)
            
            sampled_sequences = sample_rna_sequences(input_dir, k, family_name)
            temp_non_aligned_output_file = os.path.join(aln_dir, f"temp_{base_filename}")
            create_msa_fasta(sampled_sequences, temp_non_aligned_output_file)
            align_with_mafft(temp_non_aligned_output_file, aligned_output_file)
            os.remove(temp_non_aligned_output_file)  # Remove the temporary file
            format_fasta_single_line(aligned_output_file)
            print(f"Aligned MSA file created: {aligned_output_file}")

        elif align == 2:
            non_aligned_output_file = os.path.join(no_aln_dir, base_filename)
            aligned_filename = f"{family_name}.{dataset_name}.k_{k}_{i}.aln.fasta"
            aligned_output_file = os.path.join(aln_dir, aligned_filename)
            
            sampled_sequences = sample_rna_sequences(input_dir, k, family_name)
            create_msa_fasta(sampled_sequences, non_aligned_output_file)
            print(f"Non-aligned MSA file created: {non_aligned_output_file}")
            
            align_with_mafft(non_aligned_output_file, aligned_output_file)
            format_fasta_single_line(aligned_output_file)
            print(f"Aligned MSA file created: {aligned_output_file}")

def main():
    parser = argparse.ArgumentParser(description="Generate multiple MSAs from RNA sequences.")
    parser.add_argument('-i', '--input', type=str, default='./gold-database/RNAStrAlign-master/telomerase_database/',
                        help='Path to the directory containing .seq files or JSON file with arguments')
    parser.add_argument('-o', '--output_dir', type=str, default='./output/', help='Path to the directory where output files will be saved')
    parser.add_argument('-k', type=int, default=30, help='Number of sequences per MSA')
    parser.add_argument('-c', '--count', type=int, default=10, help='Number of MSAs to generate')
    parser.add_argument('-ds', '--dataset_name', type=str, default='RNAStrAlign', help='Name of the dataset')
    parser.add_argument('-f', '--family_name', type=str, default='telomerase', help='Name of the family')
    parser.add_argument('-sf', '--subfamily_name', type=str, default='telomerase', help='Name of the subfamily to include as a prefix in sequence IDs')
    parser.add_argument('-a', '--align', type=int, default=0, choices=[0, 1, 2], help='Alignment option: 0 (no alignment), 1 (do alignment), 2 (both versions)')
    
    args = parser.parse_args()
    
    if os.path.isdir(args.input):
        generate_msas(args.input, args.output_dir, args.k, args.count, args.dataset_name, args.family_name, args.subfamily_name, args.align)
    elif os.path.isfile(args.input) and args.input.endswith('.json'):
        with open(args.input, 'r') as f:
            args_list = json.load(f)
        for arg_set in args_list:
            counts = arg_set.pop('count', [10])
            output_dirs = arg_set.pop('output_dir', [args.output_dir] * len(counts))
            for count, output_dir in zip(counts, output_dirs):
                generate_msas(output_dir=output_dir, count=count, **arg_set)
    else:
        raise ValueError("Input must be a directory or a JSON file")

if __name__ == "__main__":
    main()
