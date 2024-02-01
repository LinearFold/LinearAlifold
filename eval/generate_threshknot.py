from collections import defaultdict
from covid.evaluation import get_bpp_matrix, pairs_to_struc, parse_secondary_structure, evaluate
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--threshold", type=float, default=0.3)
parser.add_argument("-i", "--bpp_file_path", type=str, required=True)
parser.add_argument("-o", "--output_file_path", type=str, required=True)
args = parser.parse_args()

def get_threshknot_pairs(bpp_file_path, threshold=0.3):
    bpp_matrix, seq_length = get_bpp_matrix(bpp_file_path, threshold)
    best_prob = defaultdict(float)

    for key, value in bpp_matrix.items():
        i, j = key
        prob = value

        best_prob[i] = max(best_prob[i], prob)
        best_prob[j] = max(best_prob[j], prob)

    pairs = []
    visited = set()
    for key, value in bpp_matrix.items():
        i, j = key
        prob = value

        if prob == best_prob[i] and prob == best_prob[j]:
            if i in visited or j in visited:
                continue

            visited.add(i)
            visited.add(j)
            pairs.append((i, j))

    return pairs, seq_length


gold_struc = huston_Arch3 = "......(((((.(((((....)))))..)))))...........(((((.....)))))[[[[[[[[.[[.[[[[[[[[[...((((((((.((.((((.(((.....))).)))))).))))))))......((((.....))))...(((((((((((..(((((...(((.(((((((((((..((((((.(((((......)))))..))))))......)))(((((((.((......)))))))))(((....)))))))))))))).))))).))))...))))))).......((((((...........((((((...))))))....)))))).....(((((.(((((((((((((.....)))).))))..))))).))))).........................((((((.((..((((......))))..))...)))))).......((((((..((.(((((((((((..(((...((......))...))).)))))))))))))))))))..............((((((((((...........))))))))).).......((((((((((......(((.........(((((((((............(((..((.(((((((((....((.((...((......))..))))...))))).))))...))....)))...............))))))))).........))).......))))))......))))..]]]].]].]]]]].]]]]]]]]..................................."

if __name__ == '__main__':
    pairs, seq_length = get_threshknot_pairs(args.bpp_file_path, args.threshold)
    pseudo_pairs = []
    pseudo_pairs2 = []
    pseudo_pairs3 = []


    # check for pseudoknots 1
    for i in range(len(pairs)):
        for j in range(i + 1, len(pairs)):
            if pairs[i][0] < pairs[j][0] < pairs[i][1] < pairs[j][1]:
                pseudo_pairs.append((pairs[j][0], pairs[j][1]))

    pseudo_pairs = set(pseudo_pairs)
    pairs = list(set(pairs) - pseudo_pairs)
    pseudo_pairs = list(pseudo_pairs)

    # check for pseudoknots 2
    for i in range(len(pseudo_pairs)):
        for j in range(i + 1, len(pseudo_pairs)):
            if pseudo_pairs[i][0] < pseudo_pairs[j][0] < pseudo_pairs[i][1] < pseudo_pairs[j][1]:
                pseudo_pairs2.append((pseudo_pairs[j][0], pseudo_pairs[j][1]))

    pseudo_pairs2 = set(pseudo_pairs2)
    pseudo_pairs = list(set(pseudo_pairs) - pseudo_pairs2)
    pseudo_pairs2 = list(pseudo_pairs2)

    # check for pseudoknots 3
    for i in range(len(pseudo_pairs2)):
        for j in range(i + 1, len(pseudo_pairs2)):
            if pseudo_pairs2[i][0] < pseudo_pairs2[j][0] < pseudo_pairs2[i][1] < pseudo_pairs2[j][1]:
                pseudo_pairs3.append((pseudo_pairs2[j][0], pseudo_pairs2[j][1]))

    pseudo_pairs3 = set(pseudo_pairs3)
    pseudo_pairs2 = list(set(pseudo_pairs2) - pseudo_pairs3)
    pseudo_pairs3 = list(pseudo_pairs3)


    print(len(pairs), len(pseudo_pairs), len(pseudo_pairs2), len(pseudo_pairs3))
            

    struc = ['.' for _ in range(seq_length)]

    for i, j in pairs:
        struc[i] = "("
        struc[j] = ")"

    for i, j in pseudo_pairs:
        struc[i] = "["
        struc[j] = "]"

    for i, j in pseudo_pairs2:
        struc[i] = "{"
        struc[j] = "}"

    for i, j in pseudo_pairs3:
        struc[i] = "<"
        struc[j] = ">"

    struc = "".join(struc)
    # print(struc)
    with open(args.output_file_path, 'w') as f:
        f.write(struc)