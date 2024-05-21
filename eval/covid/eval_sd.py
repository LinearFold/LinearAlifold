import argparse
import os
import numpy as np
import sys

from collections import defaultdict

sys.path.insert(0, os.path.abspath(os.path.join(__file__, *(['..'] * 3))))
import utility as evaluation

# hybrid reference structure
#......(((((.(((((....)))))..)))))...........(((((.....)))))[[[[[[[[.[[.[[[[[[[[[...((((((((.((.((((.(((.....))).)))))).))))))))......((((.....))))...(((((((((((..(((((...(((.(((((((((((..((((((.(((((......)))))..))))))......)))(((((((.((......)))))))))(((....)))))))))))))).))))).))))...))))))).......((((((...........((((((...))))))....)))))).....(((((.(((((((((((((.....)))).))))..))))).)))))......*..........*.........((((((.((..((((......))))..))...)))))).......((((((..((.(((((((((((..(((...((......))...))).)))))))))))))))))))..............((((((((((...........))))))))).).......((((((((((......(((.........(((((((((............(((..((.(((((((((....((.((...((......))..))))...))))).))))...))....)))...............))))))))).........))).......))))))......))))..]]]].]].]]]]].]]]]]]]]...................................

VALID_PAIRS = set(list(["AU", "UA", "CG", "GC", "UG", "GU"]))
LTF_k_25_TK = "......(((((.(((((....)))))..)))))...........(((((.....)))))((((((((.((.(((((((((((.((((((((.((.((((.(((.....))).)))))).))))))))..((((((.....)))))).(.(((((((((((..(((((...(((.((((((((.....((((((.(((((......)))))..)))))).........(((((((.((......)))))))))(((....)))))))))))))).))))).))))...))))))).).....((((((.......((..((((((...))))))..)))))))).....(((((.(((((((((((((.....)))).)))..)))))).))))).........................((((((.((..((((......))))..))...)))))).......((((((..((.(((((((((((..(((...((......))...))).))))))))))))))))))).((((....)))).((((((((((.((....))..))))))))).).......((((((((((.(((..(((.(((.....(((((((((....((......((((.((.(((((((((.(..((.((.(((.......))).)))).).))))).))))...))..).)))........)).....))))))))).....))).)))..)).).)))))......))))))))))).)).))))).))))))))..................................."
huston_Arch3_ys_ziv = "......(((((.(((((....)))))..)))))...........(((((.....)))))[[[[[[[[.[[.[[[[[[[[[...((((((((.((.((((.(((.....))).)))))).))))))))......((((.....))))...(((((((((((..(((((...(((.(((((((((((..((((((.(((((......)))))..))))))......)))(((((((.((......)))))))))(((....)))))))))))))).))))).))))...))))))).......((((((...........((((((...))))))....)))))).....(((((.(((((((((((((.....)))).))))..))))).))))).........................((((((.((..((((......))))..))...)))))).......((((((..((.(((((((((((..(((...((......))...))).)))))))))))))))))))..............((((((((((...........))))))))).).......((((((((((......(((.........(((((((((............(((..((.(((((((((....((.((...((......))..))))...))))).))))...))....)))...............))))))))).........))).......))))))......))))..]]]].]].]]]]].]]]]]]]]..................................."
huston_Arch3_no_ziv = "......(((((.(((((....)))))..)))))...........(((((.....))))).((((.......))))........((((((((.((.((((.(((.....))).)))))).))))))))......................(((((((((((..(((((...(((.(((((((((((..((((((.(((((......)))))..))))))......)))(((((((.((......)))))))))(((....)))))))))))))).))))).))))...))))))).......((((((...........((((((...))))))....)))))).....(((((.(((((((((((((.....)))).)))..)))))).))))).........................((((((.((..((((......))))..))...)))))).......((((((..((.(((((((((((..(((...((......))...))).)))))))))))))))))))..............((((((((((...........))))))))).).((((..((((((((((......(((.........(((((((((............(((..((.(((((((((....((.((...((......))..))))...))))).))))...))....)))...............))))))))).........))).......))))))......))))..))))....................................................."

parser = argparse.ArgumentParser()
parser.add_argument(
    "-d", "--data-path", type=str, default="./../data/v1", help="path to the directory containing the MSAs"
)
parser.add_argument(
    "-p",
    "--pred-path",
    type=str,
    default="./outputs/LinearAliFold_MFE_EM2_multi_detail",
    help="path to the directory containing the predicted structures",
)
parser.add_argument("-g", "--gold", type=int, default="2", help="0: LTF_k25, 1: huston_arch3, 2: huston_arch3_hybrid_ziv")
parser.add_argument("--adjust-length", action="store_true", help="whether to adjust the length of the predicted structure to match the length of the aligned sequence, only when the length of the predicted structure is less than the length of the aligned sequence")
parser.add_argument("--skip-seq", action="store_true", help="whether to skip the sequence whose predicted output file is not present")
parser.add_argument("--ignore-k1", action="store_true", help="whether to ignore the single sequence alignment family")
parser.add_argument("--backsearch", action="store_true", help="whether to use backsearch to find the predicted structure")
parser.add_argument("--slip", action="store_true", help="whether to allow slip")


def get_covid_structural_distance(data_path, pred_path, gold_id, backsearch=False, adjust_length=False, skip_seq=False, ignore_k1=False, middle_restriction=False):
    if gold_id == 0:
        gold_structure = LTF_k_25_TK
    elif gold_id == 1:
        gold_structure = huston_Arch3_no_ziv
    elif gold_id == 2:
        gold_structure = huston_Arch3_ys_ziv
    else:
        raise NotImplementedError

    families = set()
    precision, sensitivity, structural_distance, f1_score = (
        defaultdict(list),
        defaultdict(list),
        defaultdict(list),
        defaultdict(list),
    )
    five_end_utr_range = set(range(0, 401))             # 0 to 400        (both end inclusive)
    three_end_utr_range = set(range(29493, 30000))      # 29493 to 29999  (both end inclusive)

    msa_files_name = sorted([f for f in os.listdir(data_path) if (f.endswith(".fasta") or f.endswith(".txt"))])
    struc_files_name = sorted([f for f in os.listdir(pred_path) if (f.endswith(".txt"))])

    if not skip_seq:
        assert len(msa_files_name) == len(
            struc_files_name
        ), "Number of MSA (input) files and structure (output) files must be equal\nMSA files: %d\nStructure files: %d" % (
            len(msa_files_name),
            len(struc_files_name),
        )
    msa_files_name.sort(key=lambda x: int(x.split("_")[1]))
    struc_files_name.sort(key=lambda x: int(x.split("_")[1]))
    for fid in range(len(msa_files_name)):
        family = "_".join(msa_files_name[fid].split("_")[:2])
        msa_file = open(os.path.join(data_path, msa_files_name[fid])).readlines()

        # if struc file is not present, skip the sequence
        tmp_file_name = msa_files_name[fid] + ".txt"
        if tmp_file_name not in struc_files_name:
            print("Skipping sequence %s" % msa_files_name[fid])
            continue
        struc_file = open(os.path.join(pred_path, tmp_file_name)).readlines()

        for i in range(len(msa_file)):
            if msa_file[i].startswith(">Wuhan") or "Wuhan" in msa_file[0]:
                aligned_sequence = msa_file[i + 1].strip()
                break

        predicted_structure = None

        # finds the index of the predicted structure in the file that contains the output of the folding algorithm
        if backsearch:
            find_struc_index = lambda pred_file: next(
                (len(pred_file) - idx - 1 for idx, line in enumerate(reversed(pred_file)) if line.strip() and line.strip().split()[0].startswith((".", "("))), None
            )
        else:
            find_struc_index = lambda pred_file: next(
                (idx for idx, line in enumerate(pred_file) if line.strip() and line.strip().split()[0].startswith((".", "("))), None
            )

        if find_struc_index(struc_file) is not None:
            predicted_structure = struc_file[find_struc_index(struc_file)].strip().split()[0]
        else:
            print("[Error] could not find predicted structure in file %s" % struc_files_name[fid])
            raise NotImplementedError
        
        if adjust_length and len(predicted_structure) < len(aligned_sequence):
            # print("[WARNING] Length of aligned sequence and predicted structure must be equal %s" % struc_files_name[fid])
            # print("Adjusting length of the predicted structure to match the length of the aligned sequence")
            predicted_structure = predicted_structure + "." * (len(aligned_sequence) - len(predicted_structure))

        assert len(aligned_sequence) == len(
            predicted_structure
        ), "Length of aligned sequence and predicted structure must be equal %s\nSequence Length: %s\nPredicted Structure Length: %s" % (struc_files_name[fid], len(aligned_sequence), len(predicted_structure))

        a2s_map = evaluation.get_alignment_to_sequence_mapping(
            aligned_sequence
        )  # maps the indices of the aligned sequence to the indices of the unaligned sequence
        seq = aligned_sequence.replace("-", "")  # unaligned sequence

        # get the structure corresponding to the unaligned sequence
        struc = ["."] * len(seq)
        # print(predicted_structure)
        for i, j, (b1, b2) in evaluation.parse_secondary_structure(predicted_structure, True)[0]:
            if aligned_sequence[i] + aligned_sequence[j] not in VALID_PAIRS:
                continue
            if seq[a2s_map[i]] + seq[a2s_map[j]] not in VALID_PAIRS:
                continue
            if (a2s_map[i] in five_end_utr_range or a2s_map[i] in three_end_utr_range) and (
                a2s_map[j] in five_end_utr_range or a2s_map[j] in three_end_utr_range
            ):
                if a2s_map[i] == max(five_end_utr_range) or a2s_map[j] == max(five_end_utr_range):
                    continue
                struc[a2s_map[i]] = b1
                struc[a2s_map[j]] = b2                
            
            if middle_restriction and ((
                a2s_map[i] in five_end_utr_range
                and a2s_map[j] not in five_end_utr_range
                and a2s_map[j] not in three_end_utr_range
            ) or (
                a2s_map[j] in three_end_utr_range
                and a2s_map[i] not in five_end_utr_range
                and a2s_map[i] not in three_end_utr_range
            )):
                struc[a2s_map[i]] = "*"
                struc[a2s_map[j]] = "*"

        struc = "".join(struc)
        struc = (
            struc[min(five_end_utr_range) : max(five_end_utr_range)]
            + ".........."
            + struc[min(three_end_utr_range) : max(three_end_utr_range)]
        )

        # print(msa_files_name[fid])
        current_precision, current_sensitivity, _, current_structural_distance = evaluation.evaluate(
            struc, gold_structure, allow_slip=False
        )

        families.add(family)
        precision[family].append(current_precision)
        sensitivity[family].append(current_sensitivity)
        f1_score[family].append((2 * current_precision * current_sensitivity / (current_precision + current_sensitivity)) * 100)
        structural_distance[family].append(current_structural_distance)

    return families, precision, sensitivity, f1_score, structural_distance

if __name__ == "__main__":
    args = parser.parse_args()
    
    families, precision, sensitivity, f1_score, structural_distance = get_covid_structural_distance(args.data_path, args.pred_path, args.gold, args.backsearch, args.adjust_length, args.skip_seq, args.ignore_k1)

    results = []
    total_f1_avg, total_sd_avg = 0, 0
    
    if 'k_1' in structural_distance:
        for _ in range(9):
            structural_distance['k_1'].append(structural_distance['k_1'][0])

    if args.ignore_k1:
        families.remove('k_1')

    for family in sorted(families, key=lambda x: int(x.split("_")[1])):    
        results.append(
        "[{:<5}]\t{:>6.2f}\t\t{:>6.2f}\t\t({:>6.2f}  {:>6.2f}  {:>6.2f})\t[{}]\t({:>6.2f}  {:>6.2f}  {:>6.2f})".format(
            family,
            np.mean(precision[family]) * 100,
            np.mean(sensitivity[family]) * 100,
            min(f1_score[family]),
            np.mean(f1_score[family]),
            max(f1_score[family]),
            "  ".join("{:>6.2f}".format(sd) for sd in structural_distance[family]),
            min(structural_distance[family]),
            np.mean(structural_distance[family]),
            max(structural_distance[family]),
        )
    )
        total_f1_avg += np.mean(f1_score[family])
        total_sd_avg += np.mean(structural_distance[family])

    total_f1_avg /= len(families)
    total_sd_avg /= len(families)

    print()     # print an empty line
    # print(round(total_f1_avg, 3), round(total_sd_avg, 3))

    header = "[{:<5}]{:>10}\t{:>12} {:>29} {:>50}".format(
        "Family", "Precision", "Sensitivity", "F1-score (min, avg, max)", "Structural Distance (min, avg, max)"
    )

    print(header)
    print("\n".join(results))

    print("Average F1-score: %0.2f" % (total_f1_avg))
    print("Average Structural Distance: %0.2f\n" % (total_sd_avg))
