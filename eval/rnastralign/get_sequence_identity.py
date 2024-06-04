import os
import sys
import numpy as np
from collections import defaultdict

sys.path.insert(0, os.path.abspath(os.path.join(__file__, *(['..'] * 2))))
import utility


def main(data_path):
    """
    Calculate and print the average sequence identity for each MSA in the given directory.
    """
    seq_files = [f for f in os.listdir(data_path) if f.endswith(".fasta") or f.endswith(".txt")]
    seq_identities = defaultdict(list)
    seq_lens = defaultdict(list)

    # Loop through each file and calculate the sequence identities
    for file_name in seq_files:
        seqs = []
        family = file_name.split(".")[0]
        with open(os.path.join(data_path, file_name), 'r') as file:
            for line in file:
                if not line.startswith(">"):
                    seqs.append(line.strip().upper())
                    seq_lens[family].append(len(line.strip().replace("-", "")))

        family = file_name.split(".")[0]
        if len(seqs) > 1:
            avg_identity = utility.calculate_msa_seq_identity(seqs)
            seq_identities[family].append(avg_identity)

    # Print formatted results
    print("Sequence Identity")
    print("{:<8}\t{:>10}\t{:>10}".format("Family", "Identity", "Length"))

    for family in sorted(seq_identities):
        avg_value = np.mean(seq_identities[family])
        avg_len = np.mean(seq_lens[family])
        print("{:<8}\t{:>10.2f}\t{:>10.2f}".format(family, avg_value, avg_len))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data-path", type=str, required=True, help="Path to the directory containing MSA files")
    args = parser.parse_args()
    
    main(args.data_path)