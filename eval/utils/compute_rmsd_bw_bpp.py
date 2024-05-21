import matplotlib.pyplot as plt
import argparse
import math

def read_bpp(file_path):
    bpp = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) == 3:
                i, j, prob = int(parts[0]), int(parts[1]), float(parts[2])
                bpp[(i, j)] = prob
    return bpp

def compute_rmsd(bpp1, bpp2):
    rmsd = 0
    deviations = []

    all_pairs = set(bpp1.keys()).union(set(bpp2.keys()))
    intersection = set(bpp1.keys()).intersection(set(bpp2.keys()))
    print("Number of pairs in both sets: ", len(intersection))
    for pair in all_pairs:
        p1 = bpp1.get(pair, 0)
        p2 = bpp2.get(pair, 0)
        rmsd += (p1 - p2) ** 2
        deviations.append(abs(p1 - p2))
    count = len(all_pairs)
    rmsd = math.sqrt(rmsd / count)

    print("Avg deviation: ", sum(deviations) / len(deviations))
    print("Max deviation: ", max(deviations))
    print("Min deviation: ", min(deviations))

    # print box plot for deviations
    plt.ylim(0, 0.002)
    plt.boxplot(deviations, whis=[10, 90])
    plt.show()

    return rmsd

def main():
    parser = argparse.ArgumentParser(description="Compute RMSD between two sets of BPPs")
    parser.add_argument("file_path1", type=str, help="Path to the first BPP file")
    parser.add_argument("file_path2", type=str, help="Path to the second BPP file")
    args = parser.parse_args()

    bpp1 = read_bpp(args.file_path1)
    bpp2 = read_bpp(args.file_path2)

    rmsd = compute_rmsd(bpp1, bpp2)
    print(f"RMSD: {rmsd:.10}")

if __name__ == "__main__":
    main()
