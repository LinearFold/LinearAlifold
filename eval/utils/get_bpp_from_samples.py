import argparse
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(__file__, *(['..'] * 2))))
from utility import parse_secondary_structure


def count_base_pairs(file_path):
    pair_counts = {}
    total_structures = 0

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('.') or line.startswith('('):
                structure = line.strip()
                pairs, _ = parse_secondary_structure(structure)
                total_structures += 1
                for i, j in pairs:
                    assert i < j, f"Invalid pair ({i}, {j})"
                    if (i, j) not in pair_counts:
                        pair_counts[(i, j)] = 0
                    pair_counts[(i, j)] += 1

    return pair_counts, total_structures

def compute_probabilities(pair_counts, total_structures):
    probabilities = {}
    for pair, count in pair_counts.items():
        probabilities[pair] = count / total_structures
    return probabilities

def main():
    parser = argparse.ArgumentParser(description="Compute Base Pairing Probabilities from RNA Secondary Structures")
    parser.add_argument("file_path", type=str, help="Path to the file containing RNA secondary structures")
    args = parser.parse_args()

    pair_counts, total_structures = count_base_pairs(args.file_path)
    probabilities = compute_probabilities(pair_counts, total_structures)

    # Print probabilities
    for pair, prob in sorted(probabilities.items()):
        assert pair[0] < pair[1], f"Invalid pair ({pair[0]}, {pair[1]})"
        print(f"{pair[0]+1} {pair[1]+1} {prob:.5e}")

if __name__ == "__main__":
    main()
