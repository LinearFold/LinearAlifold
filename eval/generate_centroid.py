from covid.evaluation import get_bpp_matrix
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--threshold", type=float, default=0.3)
parser.add_argument("-i", "--bpp_file_path", type=str, required=True)
parser.add_argument("-o", "--output_file_path", type=str, required=True)
args = parser.parse_args()

def predict(bpp_matrix, alignment_length, bpp_threshold = 0.1, gamma = 1.0, TURN = 3):
    dp_matrix = [defaultdict(float) for _ in range(alignment_length)]
    bpp = [{} for _ in range(alignment_length)]
    paired_base = [[] for _ in range(alignment_length)]
    back_pointer = [{} for _ in range(alignment_length)]

    print("Processing BPP matrix")
    # Process BPP matrix
    for (i, j), score in list(bpp_matrix.items()):
        if score > bpp_threshold:
            bpp[i][j] = score
            paired_base[i].append(j)

    for i in range(alignment_length):
        paired_base[i].sort()

    # Dynamic programming algorithm
    for l in range(TURN + 1, alignment_length):
        print(l)
        for i in range(alignment_length - l):
            j = i + l

            dp_matrix[i][j] = dp_matrix[i + 1][j]
            for k in paired_base[i]:
                if k > j:
                    break

                score_kp1_j = 0.0
                if k < j:
                    score_kp1_j = dp_matrix[k + 1][j]

                temp_score = (gamma + 1) * bpp[i].get(k, 0) - 1 + dp_matrix[i + 1][k - 1] + score_kp1_j
                if dp_matrix[i][j] < temp_score:
                    dp_matrix[i][j] = temp_score
                    back_pointer[i][j] = k

    return dp_matrix, back_pointer

if __name__ == '__main__':
    bpp_matrix, seq_length = get_bpp_matrix(args.bpp_file_path, args.threshold)
    print(seq_length)
    dp_matrix, back_pointer = predict(bpp_matrix, seq_length)
    print(dp_matrix)
    # with open(args.output_file_path, 'w') as f:
    #     f.write(struc)


