import numpy as np
from collections import defaultdict
import argparse
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(__file__, *(['..'] * 3))))
import utility as evaluation


LTF_k_25_TK = "......(((((.(((((....)))))..)))))...........(((((.....)))))((((((((.((.(((((((((((.((((((((.((.((((.(((.....))).)))))).))))))))..((((((.....)))))).(.(((((((((((..(((((...(((.((((((((.....((((((.(((((......)))))..)))))).........(((((((.((......)))))))))(((....)))))))))))))).))))).))))...))))))).).....((((((.......((..((((((...))))))..)))))))).....(((((.(((((((((((((.....)))).)))..)))))).))))).........................((((((.((..((((......))))..))...)))))).......((((((..((.(((((((((((..(((...((......))...))).))))))))))))))))))).((((....)))).((((((((((.((....))..))))))))).).......((((((((((.(((..(((.(((.....(((((((((....((......((((.((.(((((((((.(..((.((.(((.......))).)))).).))))).))))...))..).)))........)).....))))))))).....))).)))..)).).)))))......))))))))))).)).))))).))))))))..................................."
huston_Arch3_ys_ziv = "......(((((.(((((....)))))..)))))...........(((((.....)))))[[[[[[[[.[[.[[[[[[[[[...((((((((.((.((((.(((.....))).)))))).))))))))......((((.....))))...(((((((((((..(((((...(((.(((((((((((..((((((.(((((......)))))..))))))......)))(((((((.((......)))))))))(((....)))))))))))))).))))).))))...))))))).......((((((...........((((((...))))))....)))))).....(((((.(((((((((((((.....)))).))))..))))).))))).........................((((((.((..((((......))))..))...)))))).......((((((..((.(((((((((((..(((...((......))...))).)))))))))))))))))))..............((((((((((...........))))))))).).......((((((((((......(((.........(((((((((............(((..((.(((((((((....((.((...((......))..))))...))))).))))...))....)))...............))))))))).........))).......))))))......))))..]]]].]].]]]]].]]]]]]]]..................................."
huston_Arch3_no_ziv = "......(((((.(((((....)))))..)))))...........(((((.....))))).((((.......))))........((((((((.((.((((.(((.....))).)))))).))))))))......................(((((((((((..(((((...(((.(((((((((((..((((((.(((((......)))))..))))))......)))(((((((.((......)))))))))(((....)))))))))))))).))))).))))...))))))).......((((((...........((((((...))))))....)))))).....(((((.(((((((((((((.....)))).)))..)))))).))))).........................((((((.((..((((......))))..))...)))))).......((((((..((.(((((((((((..(((...((......))...))).)))))))))))))))))))..............((((((((((...........))))))))).).((((..((((((((((......(((.........(((((((((............(((..((.(((((((((....((.((...((......))..))))...))))).))))...))....)))...............))))))))).........))).......))))))......))))..))))....................................................."


parser = argparse.ArgumentParser()
parser.add_argument(
    "-d", "--data-path", type=str, default="./../data/v1", 
    help="path to the directory containing the MSAs",
)
parser.add_argument(
    "-p", "--pred-path", type=str, 
    help="path to the directory containing the predicted structures",
)
parser.add_argument("-g", "--gold", type=int, default="2", help="0: LTF_k25, 1: huston_arch3, 2: huston_arch3_hybrid_ziv")
parser.add_argument("--zero-based-idx", action="store_true", help="use zero-based index for the paired probabilities")

# global variables
five_end_utr_range = set(range(0, 401))
join_length = 10
three_end_utr_range = set(range(29493, 30000))
utr_join_index = max(five_end_utr_range) + join_length


def get_paired_probs(aligned_sequence_file_path, base_pairing_prob_file_path, k, zero_based_idx):
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
        if zero_based_idx:
            l, r = int(l), int(r)
        else:
            l, r = int(l) - 1, int(r) - 1

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


def get_covid_ensemble_defect(seqs_path, preds_path, gold_id, zero_based_idx=True):
    if gold_id == 0:
        gold_struc = LTF_k_25_TK
    elif gold_id == 1:
        gold_struc = huston_Arch3_no_ziv
    elif gold_id == 2:
        gold_struc = huston_Arch3_ys_ziv
    else:
        raise NotImplementedError

    msa_files_name = sorted([f for f in os.listdir(seqs_path) if f.endswith(".fasta") or f.endswith(".txt")])
    prob_files_name = sorted([f for f in os.listdir(preds_path) if f.endswith(".txt")])

    msa_files_name.sort(key=lambda x: int(x.split("_")[1]))
    prob_files_name.sort(key=lambda x: int(x.split("_")[1]))

    # for i in range(len(msa_files_name)):
        

    # assert len(msa_files_name) == len(
    #     prob_files_name
    # ), "Number of MSA (input) files and structure (output) files must be equal"

    gold_paired_non_adjusted, gold_unpaired_non_adjusted = evaluation.parse_secondary_structure(gold_struc)
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

    families = set()
    ensemble_defect_family = defaultdict(list)
    count_family = defaultdict(int)

    for fid in range(len(msa_files_name)):
        family = "_".join(msa_files_name[fid].split("_")[:2])
        msa_file_path = os.path.join(seqs_path, msa_files_name[fid])

        # if prob file is not present, skip the sequence
        tmp_file_name = msa_files_name[fid] + ".txt"
        if tmp_file_name not in prob_files_name:
            print("Skipping sequence %s" % msa_files_name[fid])
            continue
        prob_file_path = os.path.join(preds_path, tmp_file_name)

        # check if prob_file_path is valid
        if not os.path.exists(prob_file_path):
            print("[Warning] %s does not exist" % prob_file_path)
            continue
        
        paired_probs, seq_length = get_paired_probs(msa_file_path, prob_file_path, int(family.split("_")[1]), zero_based_idx)
        unpaired_probs = get_unpaired_probs(paired_probs, seq_length)
        ensemble_defect_family[family].append(
            round(
                compute_ensemble_defect(
                    gold_paired, gold_unpaired, len(gold_struc) - join_length, paired_probs, unpaired_probs
                ),
                2,
            )
        )
        count_family[family] += 1
        families.add(family)

    return families, count_family, ensemble_defect_family


if __name__ == "__main__":
    args = parser.parse_args()

    families, count_family, ensemble_defect_family = get_covid_ensemble_defect(args.data_path, args.pred_path, args.gold, args.zero_based_idx)

    results = []
    if "k_1" in ensemble_defect_family:
        for _ in range(9):
            ensemble_defect_family["k_1"].append(ensemble_defect_family["k_1"][0])

    for family in sorted(ensemble_defect_family.keys(), key=lambda x: int(x.split("_")[1])):
        results.append(
            "[{:<5}]\t {}\t\t{:0.2f}".format(
                family, "  ".join("{:>6.2f}".format(sd) for sd in ensemble_defect_family[family]), np.mean(ensemble_defect_family[family])
            )
        )

    print("[Family] \tEnsemble Defect\t Avg")
    print("\n".join(results))
   