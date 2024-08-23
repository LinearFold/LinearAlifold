# Example run commands:
# python3 eval_perf.py -d ./data/no_aln/ -p ./outputs/lp_v_mea --allow-slip
# python3 eval_perf.py -d ./data/aln/ -p ./outputs/laf_part_vn_ma_1.2_0.1/threshknot --cnsns --allow-slip
# python3 eval_perf.py -d ./data/aln/ -p ./outputs/raf_mfe_mea_rna/output --cnsns --allow-slip --backsearch

import argparse
import os
from collections import defaultdict
import sys
import numpy as np


sys.path.insert(0, os.path.abspath(os.path.join(__file__, *(['..'] * 2))))
import utility

def evaluate_rnastralign_performance(data_path, pred_path, ct_path = "./data/gold-database/", cnsns=False, get_msa_seq_info=False, \
                                     skip_seq=True, allow_slip=True, verbose=False, backsearch=False):
    seq_files_name = sorted([f for f in os.listdir(data_path) if (f.endswith(".fasta") or f.endswith(".txt"))])
    struc_files_name = sorted([f for f in os.listdir(pred_path) if (f.endswith(".txt"))])

    if not skip_seq:
        assert len(seq_files_name) == len(
            struc_files_name
        ), "Number of MSA (input) files and structure (output) files must be equal\nMSA files: %d\nStructure files: %d" % (
            len(seq_files_name),
            len(struc_files_name),
        )

    families = set()
    precision_scores = defaultdict(lambda: list())
    sensitivity_scores = defaultdict(lambda: list())
    f1_scores = defaultdict(lambda: list())
    structural_distances = defaultdict(lambda: list())
    seq_info = defaultdict(lambda: list())
    

    for fid in range(len(seq_files_name)):
        if len(seq_files_name) != len(struc_files_name) and seq_files_name[fid] + ".txt" not in struc_files_name:
            print("Skipping", seq_files_name[fid], "as structure file not found")
            continue
        
        if verbose:
            print("Processing", seq_files_name[fid])
        family = seq_files_name[fid].split(".")[0]
        if family == "23S":
            if verbose:
                print("Skipping 23S")
            continue

        families.add(family)
        seqs, ref_data, pred_strucs = [], [], []
        for line in open(os.path.join(data_path, seq_files_name[fid]), "r").readlines():
            if line.startswith(">"):
                ref_struc_path = os.path.join(ct_path, line.strip()[1:])
                ref_data.append(utility.parse_ct_file(ref_struc_path))
            else:
                seqs.append(line.strip().upper())

        # calculate sequence identity
        if get_msa_seq_info:
            for seq in seqs:
                assert all([x in "ACGUN-" for x in seq]), "Sequence contains non ACGU- characters: " + str(set(seq))

            msa_seq_identity = utility.calculate_msa_seq_identity(seqs)
            msa_seq_avg_len = np.mean([len(seq.replace("-", "")) for seq in seqs])
            seq_info[family].append({"identity": msa_seq_identity, "avg_len": msa_seq_avg_len})

        struc_lines = open(os.path.join(pred_path, struc_files_name[fid]), "r").readlines()
        if backsearch:
            struc_lines = struc_lines[::-1]
        for line in struc_lines:
            if line[0] in set(["(", ".", "<", "[", "{"]):
                if cnsns:
                    consns_struc = line.strip().split()[0]
                    for seq in seqs:
                        pred_strucs.append(utility.map_consns_struc_to_aln_seq(consns_struc, seq)[0])
                    break
                else:
                    pred_strucs.append(line.strip())

        assert len(seqs) == len(pred_strucs), "Number of sequences and structures must be equal " + str(len(seqs)) + " " + str(len(pred_strucs))
        assert len(ref_data) == len(pred_strucs), "Number of reference structures and structures must be equal " + str(len(ref_data)) + " " + str(len(pred_strucs))

        for pred_struc, (ref_paired, ref_unpaired) in zip(pred_strucs, ref_data):
            ref_len = 2 * len(ref_paired) + len(ref_unpaired)
            assert len(pred_struc) == ref_len, "Length of predicted structure and reference structure must be equal. " + str(len(pred_struc)) + " " + str(ref_len)
            precision, sensitivity, f1, structural_distance = utility.evaluate(pred_struc, gold_paired_pos_tuple=ref_paired, gold_unpaired_pos=ref_unpaired, allow_slip=allow_slip)
            precision_scores[family].append(precision)
            sensitivity_scores[family].append(sensitivity)
            f1_scores[family].append(f1)
            structural_distances[family].append(structural_distance)

    return families, precision_scores, sensitivity_scores, f1_scores, structural_distances, seq_info


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data-path", type=str, default="./data/no_aln", help="path to data folder")
    parser.add_argument("-p", "--pred-path", type=str, default="./outputs", help="path to prediction folder")
    parser.add_argument("-ct", "--ct-path", type=str, default="./data/gold-database/", help="path to database folder")
    parser.add_argument("--skip-seq", action="store_true", help="skip sequence file if structure file is not found")
    parser.add_argument("--verbose", action="store_true", help="print verbose output")
    parser.add_argument("--cnsns", action="store_true", help="predicted structure is a single consensus structure")
    parser.add_argument("--allow-slip", action="store_true", help="allow slip in base pair matching")
    parser.add_argument("--seq-identity", action="store_true", help="calculate the sequence identity")
    parser.add_argument("--backsearch", action="store_true", help="search for the structure in reverse order")
    args = parser.parse_args()

    train_set = set(["tRNA", "5S", "tmRNA", "group"])

    families, precision_scores, sensitivity_scores, f1_scores, structural_distances, seq_info = \
        evaluate_rnastralign_performance(args.data_path, args.pred_path, args.ct_path, skip_seq=args.skip_seq, backsearch=args.backsearch,
                                         verbose=args.verbose, cnsns=args.cnsns, allow_slip=args.allow_slip, get_msa_seq_info=args.seq_identity)

    total_f1_avg = 0
    total_sd_avg = 0
    
    results = []
    for family in sorted(families):
        results.append(
        "{:<8}\t{:>6.2f}\t\t{:>6.2f}\t\t({:>4.2f}  {:>6.2f}  {:>6.2f})\t\t({:>6.2f}  {:>6.2f}  {:>6.2f})".format(
            family,
            np.mean(precision_scores[family]) * 100,
            np.mean(sensitivity_scores[family]) * 100,
            min(f1_scores[family]) * 100,
            np.mean(f1_scores[family]) * 100,
            max(f1_scores[family]) * 100,
            min(structural_distances[family]),
            np.mean(structural_distances[family]),
            max(structural_distances[family]),
        )
    )

        total_f1_avg += np.mean(f1_scores[family])
        total_sd_avg += np.mean(structural_distances[family])

    total_f1_avg /= len(families)
    total_sd_avg /= len(families)

    print()     # print an empty line

    header = "{:<8}\t{:>10}\t{:>10}\t{:>10}\t{:>20}".format(
        "Family", "Precision", "Sensitivity", "F1-score (min, avg, max)", "Structural Distance (min, avg, max)"
    )

    print(header)
    print("\n".join(results))

    print("Average F1-score: %0.2f" % (total_f1_avg * 100))
    print("Average Structural Distance: %0.2f\n" % (total_sd_avg))

    # show sequence identity if requested
    if args.seq_identity:
        average_value = {}
        for family in families:
            average_value[family] = np.mean(seq_info[family])
        print("Sequence Identity")
        print("{:<8}\t{:>10}".format("Family", "Identity"))
        for family in sorted(families):
            print("{:<8}\t{:>10.2f}".format(family, average_value[family]))
        print()
