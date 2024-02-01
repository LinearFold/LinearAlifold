# run command example: python3 ./run_eval.py -p  ./../../tmp/part_rnastraln_bl_md/mea/ -m 2
import argparse
import os
from collections import defaultdict

import evaluation

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--data-path", type=str, default="./data", help="path to data folder")
parser.add_argument("-p", "--pred-path", type=str, default="./outputs", help="path to prediction folder")
parser.add_argument("-ct", "--ct-path", type=str, default="./database", help="path to database folder")

parser.add_argument(
    "-o",
    "--option",
    type=int,
    default="0",
    help="0: LinAliFold, 1: LinearAliFold_MFE (EM 1), 2: LinearAliFold_MFE (EM 2, multi approx), 3: LinearAliFold_MFE (EM 2, multi detail",
)
parser.add_argument("-m", "--mode", type=int, default="0", help="0: train, 1: test, 2: overall")


pairables = ["AU", "UA", "CG", "GC", "UG", "GU"]

train_set = set(["tRNA", "5S", "tmRNA", "group"])


if __name__ == "__main__":
    args = parser.parse_args()
    data_path = args.data_path
    pred_path = args.pred_path

    files = sorted(os.listdir(data_path))

    test_seqs = defaultdict(list)
    for file in files:
        if ".fasta" not in file or "txt" in file or "bpp" in file:
            continue

        temp = open(os.path.join(data_path, file)).readlines()

        for index, line in enumerate(temp):
            line = line.strip()
            if line[0] == ">":
                test_seqs[file].append((line, temp[index + 1].strip().upper()))

    CNT = 0
    tot_P_slip, tot_R_slip, tot_P_non_slip, tot_R_non_slip = (
        defaultdict(float),
        defaultdict(float),
        defaultdict(float),
        defaultdict(float),
    )
    tot_F_slip, tot_F_non_slip = defaultdict(float), defaultdict(float)
    tot_cnt = defaultdict(int)
    all_families = set()

    tot_score, tot_FE, tot_covar = defaultdict(float), defaultdict(float), defaultdict(float)
    tot_pred_cnt = defaultdict(int)

    for file, elem in sorted(test_seqs.items(), key=lambda x: x[0]):

 
        
        pred_file = os.path.join(pred_path, file + ".txt")

        
        family = file.split("_")[0].split(".")[0]

        if args.mode == 0:
            if family not in train_set:
                continue
        elif args.mode == 1:
            if family in train_set:
                continue

        all_families.add(family)

        temp = open(pred_file).readlines()

        # if temp[0][0] == 'W':
        #     print(1)
        #     temp = temp[::-1]

        # go to line that starts with '.' or '('
        # i = 0
        # while i < len(temp) and temp[i][0] not in ".(":
        #     i += 1

        i = len(temp) - 1
        while i >= 0 and temp[i][0] not in ".(":
            i -= 1

        if i == len(temp):
            print("Error in file: {}".format(file))
            continue
        # print(family)

        pred_stru = temp[i].strip().split()[0]
        tot_score[family] += 0
        tot_FE[family] += 0
        tot_covar[family] += 0
        # scores = temp[0].strip().split()[1:]  # .split(' ')
        # tot_score[family] += float(scores[0])
        # tot_FE[family] += float(scores[3])
        # tot_covar[family] += float(scores[-1])

        # print(len(pred_stru), len(elem[0][1]))
        # if len(pred_stru) < len(elem[0][1]):
        #     pred_stru += "." * (len(elem[0][1]) - len(pred_stru))
            # print("Warning: {} is shorter than {}".format(pred_file, file))
        # exit()

        tot_pred_cnt[family] += 1

        skip = False

        for ct_path, seq in elem:
            ct_path = os.path.join(args.ct_path, ct_path[3:])
            temp = open(ct_path).readlines()

            for nn in seq:
                if nn not in "AUCG-":
                    skip = True

            if "23s" in ct_path:
                if temp[4][0] != "1":
                    raise NotImplementedError
                nucs = ""
                gold_pairs = set()
                for index, line in enumerate(temp[4:]):
                    l, nuc, r = line.strip().split()
                    l = int(l)
                    r = int(r)
                    nucs += nuc
                    if l >= r:
                        continue
                    gold_pairs.add((l, r))

                if seq.replace("-", "") != nucs:
                    raise NotImplementedError
                pred_pairs = evaluation.pairs(pred_stru)
                pred_stru_dot_bracket = ["."] * len(seq)
                for l, r in pred_pairs:
                    if seq[l - 1] == "-" or seq[r - 1] == "-" or seq[l - 1] + seq[r - 1] not in pairables:
                        # if seq[l-1] == "-":
                        #     pred_stru_dot_bracket[l-1] = '-'
                        # if seq[r-1] == '-':
                        #     pred_stru_dot_bracket[r-1] = '-'
                        continue

                    pred_stru_dot_bracket[l - 1] = "("
                    pred_stru_dot_bracket[r - 1] = ")"

                for ii, nuc in enumerate(seq):
                    if nuc == "-":
                        pred_stru_dot_bracket[ii] = "-"

                pred_stru_dot_bracket = ("".join(pred_stru_dot_bracket)).replace("-", "")

                if len(pred_stru_dot_bracket) != len(nucs):
                    print(pred_stru_dot_bracket)
                    print(seq)
                    print(nucs)
                    raise NotImplementedError

                pred_clean_pairs = evaluation.pairs(pred_stru_dot_bracket)

                P_slip, R_slip, P_non_slip, R_non_slip = evaluation.eval(gold_pairs, pred_clean_pairs)

            else:
                length = int(temp[0].strip().split()[0])

                if length != len(temp) - 1:
                    raise NotImplementedError
                nucs = ""
                gold_pairs = set()
                for index, line in enumerate(temp[1:]):
                    l, nuc, _, _, r, _ = line.strip().split()
                    l = int(l)
                    r = int(r)
                    nucs += nuc
                    if l >= r:
                        continue
                    gold_pairs.add((l, r))

                if seq.replace("-", "") != nucs:
                    raise NotImplementedError
                pred_pairs = evaluation.pairs(pred_stru)
                pred_stru_dot_bracket = ["."] * len(seq)
                for l, r in pred_pairs:
                    if seq[l - 1] == "-" or seq[r - 1] == "-" or seq[l - 1] + seq[r - 1] not in pairables:
                        # if seq[l-1] == "-":
                        #     pred_stru_dot_bracket[l-1] = '-'
                        # if seq[r-1] == '-':
                        #     pred_stru_dot_bracket[r-1] = '-'
                        continue

                    pred_stru_dot_bracket[l - 1] = "("
                    pred_stru_dot_bracket[r - 1] = ")"

                for ii, nuc in enumerate(seq):
                    if nuc == "-":
                        pred_stru_dot_bracket[ii] = "-"

                pred_stru_dot_bracket = ("".join(pred_stru_dot_bracket)).replace("-", "")

                if len(pred_stru_dot_bracket) != len(nucs):
                    print(pred_stru_dot_bracket)
                    print(seq)
                    print(nucs)
                    raise NotImplementedError

                pred_clean_pairs = evaluation.pairs(pred_stru_dot_bracket)

                P_slip, R_slip, P_non_slip, R_non_slip = evaluation.eval(gold_pairs, pred_clean_pairs)

            F_slip = 2 * P_slip * R_slip / (P_slip + R_slip + 1e-10)
            F_non_slip = 2 * P_non_slip * R_non_slip / (P_non_slip + R_non_slip + 1e-10)
            tot_P_slip[family] += P_slip
            tot_R_slip[family] += R_slip
            tot_P_non_slip[family] += P_non_slip
            tot_R_non_slip[family] += R_non_slip
            tot_F_slip[family] += F_slip
            tot_F_non_slip[family] += F_non_slip
            tot_cnt[family] += 1

        CNT += 1
        if skip:
            print(file)

    all_families = sorted(all_families)

    results_slip = []
    results_non_slip = []
    for family in all_families:
        P_slip = tot_P_slip[family] * 100 / tot_cnt[family]
        R_slip = tot_R_slip[family] * 100 / tot_cnt[family]
        F_slip = tot_F_slip[family] * 100 / tot_cnt[family]
        P_non_slip = tot_P_non_slip[family] * 100 / tot_cnt[family]
        R_non_slip = tot_R_non_slip[family] * 100 / tot_cnt[family]
        F_non_slip = tot_F_non_slip[family] * 100 / tot_cnt[family]

        # avrg_score = tot_score[family] / tot_pred_cnt[family]
        # avrg_FE = tot_FE[family] / tot_pred_cnt[family]
        # avrg_covar = tot_covar[family] / tot_pred_cnt[family]

        results_slip.append("{:<10}\t\t{:0.2f}\t{:0.2f}\t\t{:0.2f}".format(family, P_slip, R_slip, F_slip))
        results_non_slip.append(
            "{:<10}\t\t{:0.2f}\t{:0.2f}\t\t{:0.2f}".format(family, P_non_slip, R_non_slip, F_non_slip)
        )

        # print('[Family]: {}\nP_slip = {:0.2f}, R_slip = {:0.2f}, F_slip = {:0.2f}\nP_noslip = {:0.2f}, R_noslip = {:0.2f}, F_noslip = {:0.2f}'.format(
        #     family, P_slip, R_slip, F_slip, P_non_slip, R_non_slip, F_non_slip))
        # print("avrg_score = {:.2f} | avrg_FE = {:.2f} | avrg_covar = {:.2f}".format(avrg_score, avrg_FE, avrg_covar))

    # add average for all families
    P_slip = sum(tot_P_slip.values()) * 100 / sum(tot_cnt.values())
    R_slip = sum(tot_R_slip.values()) * 100 / sum(tot_cnt.values())
    F_slip = sum(tot_F_slip.values()) * 100 / sum(tot_cnt.values())
    P_non_slip = sum(tot_P_non_slip.values()) * 100 / sum(tot_cnt.values())
    R_non_slip = sum(tot_R_non_slip.values()) * 100 / sum(tot_cnt.values())
    F_non_slip = sum(tot_F_non_slip.values()) * 100 / sum(tot_cnt.values())
    results_slip.append("{:<10}\t\t{:0.2f}\t{:0.2f}\t\t{:0.2f}".format("Average", P_slip, R_slip, F_slip))
    results_non_slip.append(
        "{:<10}\t\t{:0.2f}\t{:0.2f}\t\t{:0.2f}".format("Average", P_non_slip, R_non_slip, F_non_slip)
    )

    print("------------------ Slip mode ------------------")
    print("Family\t\tPrecision\tSensitivity\t\tF1-score")
    print("\n".join(results_slip))
    # print("\n")
    # print("---------------- Non-slip mode ----------------")
    # print("Family\t\tPrecision\tSensitivity\t\tF1-score")
    # print("\n".join(results_non_slip))

    # print('[Family] Precision\tSensitivity\tF1-score\tStructural Distance')
    # print('\n'.join(results))
