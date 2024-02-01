import argparse
import evaluation
import os
from collections import defaultdict
import numpy as np

parser = argparse.ArgumentParser()
# add argument for input directory and output directory
parser.add_argument('-i', '--input', type=str, default='./sequences', help='input directory')
parser.add_argument('-o', '--output', type=str, default='./structures', help='output directory')


if __name__ == '__main__':
    args = parser.parse_args()

    input_files = os.listdir(args.input)
    output_files = os.listdir(args.output)
    assert len(input_files) == len(output_files), 'number of input files and output files do not match'
    

    input_files.sort()
    output_files.sort()

    precision, sensitivity, structural_distance, f1_score, count = (
        defaultdict(float),
        defaultdict(float),
        defaultdict(list),
        defaultdict(list),
        defaultdict(float),
    )

    i = 0
    # open each file and compare the sequence
    for input_file, output_file in zip(input_files, output_files):
        
        family = input_file.split('_')[0]
        with open(os.path.join(args.input, input_file), 'r') as f:
            for line in f.readlines():
                if line.startswith('.') or line.startswith('('):
                    gold_struc = line.strip()
                    break
            
        with open(os.path.join(args.output, output_file), 'r') as f:
            lines = f.readlines()
            for line in reversed(lines):
                if line.startswith('.') or line.startswith('('):
                    output_struc = line.strip().split()[0]
                    break

        if len(output_struc) < len(gold_struc):
            output_struc += '.' * (len(gold_struc) - len(output_struc))
    
        valid, prec, sens, f1, dist = evaluation.evaluate(output_struc, gold_struc, allow_slip=True)
        
        if valid:
            precision[family] += prec * 100
            sensitivity[family] += sens * 100
            f1_score[family].append(f1 * 100)
            structural_distance[family].append(dist)
            count[family] += 1
        else:
            print(input_file, output_struc, gold_struc)

    results = []
    total_f1_avg, total_sd_avg = 0, 0


    
    p_avg, s_avg, f_avg, sd_avg = 0, 0, 0, 0
    for family in sorted(precision.keys()):
        precision[family] = precision[family] / count[family]
        sensitivity[family] = sensitivity[family] / count[family]

        p_avg += precision[family]
        s_avg += sensitivity[family]
        f_avg += np.mean(f1_score[family])
        sd_avg += np.mean(structural_distance[family])

    p_avg = p_avg / len(precision)
    s_avg = s_avg / len(sensitivity)
    f_avg = f_avg / len(f1_score)
    sd_avg = sd_avg / len(structural_distance)
        

    for family in sorted(precision.keys()):
        # precision[family] = precision[family] / count[family]
        # sensitivity[family] = sensitivity[family] / count[family]

        
        results.append(
            "[{:<10}]\t {:>6.2f}\t\t{:>6.2f}\t\t{:>6.2f}  {:>6.2f}  {:>6.2f}\t\t{:>6.2f}   {:>6.2f}   {:>6.2f}".format(
                family,
                precision[family],
                sensitivity[family],
                min(f1_score[family]),
                np.mean(f1_score[family]),
                max(f1_score[family]),

                min(structural_distance[family]),
                np.mean(structural_distance[family]),
                max(structural_distance[family]),
            )
        )


        total_f1_avg += np.mean(f1_score[family])
        total_sd_avg += np.mean(structural_distance[family])

    results.append(
    "{:<10}\t {:>6.2f}\t\t{:>6.2f}\t\t{:>6.2f}  {:>6.2f}  {:>6.2f}\t\t{:>6.2f}   {:>6.2f}   {:>6.2f}".format(
        "Average",
        p_avg,
        s_avg,
        float("nan"),
        f_avg,
        float("nan"),
        float("nan"),
        sd_avg,
        float("nan"),
    ))

    total_f1_avg /= len(precision)
    total_sd_avg /= len(precision)

    print("{:0.2f} {:0.2f}".format(round(total_f1_avg, 2), round(total_sd_avg, 2)))

    print("[Family]\tPrecision\tSensitivity\tF1-score (min, avg, max)\tStructural Distance (min, avg, max)")
    print("\n".join(results))
    # print("Average F1-score: {:0.2f}".format(total_f1_avg))
    # print("Average Structural Distance: {:0.2f}".format(total_sd_avg))





    
    
