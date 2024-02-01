import evaluation
import numpy as np
from collections import defaultdict
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument(
    "-d", "--data-path", type=str, default="./data/v1/", help="path to the directory containing the MSAs",
)
parser.add_argument("--no-idx-offset", action="store_true", help="whether to not offset the index by 1")
parser.add_argument(
    "-p",
    "--pred-path",
    type=str,
    default="./outputs/LinearAliFold_MFE_EM2_multi_detail",
    help="path to the directory containing the predicted structures",
)
parser.add_argument("-g", "--gold", type=int, default="1", help="0: LTF_k_25_TK, 1: huston_Arch3")

# global variables
huston_Arch3 = "......(((((.(((((....)))))..)))))...........(((((.....)))))[[[[[[[[.[[.[[[[[[[[[...((((((((.((.((((.(((.....))).)))))).))))))))......((((.....))))...(((((((((((..(((((...(((.(((((((((((..((((((.(((((......)))))..))))))......)))(((((((.((......)))))))))(((....)))))))))))))).))))).))))...))))))).......((((((...........((((((...))))))....)))))).....(((((.(((((((((((((.....)))).))))..))))).))))).........................((((((.((..((((......))))..))...)))))).......((((((..((.(((((((((((..(((...((......))...))).)))))))))))))))))))..............((((((((((...........))))))))).).......((((((((((......(((.........(((((((((............(((..((.(((((((((....((.((...((......))..))))...))))).))))...))....)))...............))))))))).........))).......))))))......))))..]]]].]].]]]]].]]]]]]]]..................................."
five_end_utr_range = set(range(0, 401))
join_length = 10
three_end_utr_range = set(range(29493, 30000))
utr_join_index = max(five_end_utr_range) + join_length


def get_paired_probs(aligned_sequence_file_path, base_pairing_prob_file_path, k):
    # read file
    with open(aligned_sequence_file_path, "r") as f:
        aligned_sequence = f.readlines()
    with open(base_pairing_prob_file_path, "r") as f:
        pairing_prob_lines = f.readlines()

    # if ">NC_045512.2_Wuhan_seafood_market_pneumonia_virus_isolate_Wuhan-Hu-1__complete_genome" not in aligned_sequence[0].strip():
    #     raise NotImplementedError

    aligned_sequence = aligned_sequence[1].strip()
    ungapped_sequence = aligned_sequence.replace("-", "")

    a2s_index = evaluation.get_alignment_to_sequence_mapping(aligned_sequence)
    paired_probs = defaultdict(lambda: defaultdict(float))

    for line in pairing_prob_lines[k:]:
        if line == "\n" or len(line.split()) != 3 or line.startswith(".") or line.startswith("("):
            # print(line)
            continue

        l, r, prob = line.split()

        if not args.no_idx_offset:
            l, r = int(l) - 1, int(r) - 1
        else:
            l, r = int(l), int(r)

        if l >= len(aligned_sequence):
            print(l, r, prob)
            print(aligned_sequence_file_path, len(aligned_sequence))

        # check if l, r actually pair in the reference sequence
        if aligned_sequence[l] + aligned_sequence[r] not in evaluation.VALID_PAIRS:
            continue
        
        l, r = a2s_index[l], a2s_index[r]

        if ungapped_sequence[l] + ungapped_sequence[r] not in evaluation.VALID_PAIRS:
            continue
        if (l in five_end_utr_range or l in three_end_utr_range) and (
            r in five_end_utr_range or r in three_end_utr_range
        ):
            if l == max(five_end_utr_range) or r == max(five_end_utr_range):
                continue

            paired_probs[l][r] = float(prob)

    return paired_probs, len(ungapped_sequence)


def get_unpaired_probs(paired_probs, seq_length):
    unpaired_probs = {}

    for i in range(max(five_end_utr_range)):
        unpaired_probs[i] = 1.00
    for i in range(min(three_end_utr_range), seq_length):
        unpaired_probs[i] = 1.00

    for l in paired_probs.keys():
        for r, p in paired_probs[l].items():
            unpaired_probs[l] -= p
            unpaired_probs[r] -= p

    return unpaired_probs


def compute_ensemble_defect(gold_paired, gold_unpaired, gold_length, paired_probs, unpaired_probs):
    value = gold_length
    for pair in gold_paired:
        value -= 2 * paired_probs[pair[0]][pair[1]]
    for i in gold_unpaired:
        value -= unpaired_probs[i]
    return value


if __name__ == "__main__":
    args = parser.parse_args()
    msa_files_name = sorted([f for f in os.listdir(args.data_path) if f.endswith(".fasta")])
    prob_files_name = sorted([f for f in os.listdir(args.pred_path) if f.endswith(".txt")])

    # assert len(msa_files_name) == len(
    #     prob_files_name
    # ), "Number of MSA (input) files and structure (output) files must be equal"

    gold_paired_non_adjusted, gold_unpaired_non_adjusted = evaluation.parse_secondary_structure(huston_Arch3)
    gold_paired, gold_unpaired = [], []
    for l, r in gold_paired_non_adjusted:
        if l >= utr_join_index:
            l = min(three_end_utr_range) + l - utr_join_index
        if r >= utr_join_index:
            r = min(three_end_utr_range) + r - utr_join_index
        gold_paired.append((l, r))
    gold_paired.sort()

    for j in gold_unpaired_non_adjusted:
        if j >= max(five_end_utr_range) and j < utr_join_index:
            continue
        if j >= utr_join_index:
            j = min(three_end_utr_range) + j - utr_join_index
        gold_unpaired.append(j)

    ensemble_defect_family = defaultdict(list)
    count_family = defaultdict(int)

    for fid in range(len(msa_files_name)):

        family = "_".join(msa_files_name[fid].split("_")[:2])
        msa_file_path = os.path.join(args.data_path, msa_files_name[fid])
        prob_file_path = os.path.join(args.pred_path, prob_files_name[fid])

        # check if prob_file_path is valid
        if not os.path.exists(prob_file_path):
            print("[Warning] %s does not exist" % prob_file_path)
            continue
        
        paired_probs, seq_length = get_paired_probs(msa_file_path, prob_file_path, int(family.split("_")[1]))
        unpaired_probs = get_unpaired_probs(paired_probs, seq_length)
        ensemble_defect_family[family].append(
            round(
                compute_ensemble_defect(
                    gold_paired, gold_unpaired, len(huston_Arch3) - join_length, paired_probs, unpaired_probs
                ),
                2,
            )
        )
        count_family[family] += 1

    results = []
    for _ in range(9):
        ensemble_defect_family["k_1"].append(ensemble_defect_family["k_1"][0])
    for family in sorted(ensemble_defect_family.keys(), key=lambda x: int(x.split("_")[1])):
        results.append(
            "[{:<5}]\t {}\t\t{:0.2f}".format(
                family, ensemble_defect_family[family], np.mean(ensemble_defect_family[family])
            )
        )

    print("[Family] \tEnsemble Defect\t Avg")
    print("\n".join(results))
