from hopcroftkarp import HopcroftKarp
from collections import defaultdict

VALID_PAIRS = set(["AU", "UA", "CG", "GC", "GU", "UG"])

def parse_secondary_structure(struc, get_pair_info = False):
    """
    Parse an RNA secondary structure string and return the indices of paired and unpaired nucleotides.

    Args:
        struc (str): A string representing RNA secondary structure, where '(' and ')' denote paired nucleotides
                     and '.' denotes unpaired nucleotides.

    Returns:
        tuple: A tuple containing two lists:
            - A list of tuples representing the indices of paired nucleotides in the structure string.
            - A list of indices representing the indices of unpaired nucleotides in the structure string.

    Example:
        >>> parse_secondary_structure('((..((...)).))')
        ([(0, 11), (4, 9), (5, 8)], [2, 3, 6, 7, 12, 13])

    If the input string contains unbalanced parentheses, the function returns None and prints an error message.
    """
    stack1, stack2, stack3, stack4 = [], [], [], []
    paired_indices = []  # list of tuples: [(i1, j1), (i2, j2), ...]
    unpaired_indices = []  # list of indices: [i1, i2, ...]

    try:
        for i, x in enumerate(struc):
            if x == "(":
                stack1.append(i)
            elif x == "[":
                stack2.append(i)
            elif x == "{":
                stack3.append(i)
            elif x == "<":
                stack4.append(i)
            elif x == ")":
                u, v = stack1.pop(), i
                assert struc[u] == '(' and struc[v] == ')', "Unbalanced parenthesis in structure string"
                paired_indices.append((u, v, ('(', ')')) if get_pair_info else (u, i))
            elif x == "]":
                # paired_indices.append((stack2.pop(), i, ('[', ']')) if get_pair_info else (stack2.pop(), i))
                u, v = stack2.pop(), i
                assert struc[u] == '[' and struc[v] == ']', "Unbalanced parenthesis in structure string"
                paired_indices.append((u, v, ('[', ']')) if get_pair_info else (u, v))
            elif x == "}":
                # paired_indices.append((stack3.pop(), i, ('{', '}')) if get_pair_info else (stack3.pop(), i))
                u, v = stack3.pop(), i
                assert struc[u] == '{' and struc[v] == '}', "Unbalanced parenthesis in structure string"
                paired_indices.append((u, v, ('{', '}')) if get_pair_info else (u, v))
            elif x == ">":
                # paired_indices.append((stack4.pop(), i, ('<', '>')) if get_pair_info else (stack4.pop(), i))
                u, v = stack4.pop(), i
                assert struc[u] == '<' and struc[v] == '>', "Unbalanced parenthesis in structure string"
                paired_indices.append((u, v, ('<', '>')) if get_pair_info else (u, v))
            elif x == ".":
                unpaired_indices.append(i)
    except Exception as _:
        print("[Error] Unbalanced parenthesis in structure string")
        return None

    if stack1 or stack2 or stack3 or stack4:
        print("[Error] Unbalanced parenthesis in structure string")

    return paired_indices, unpaired_indices


def evaluate(pred_struc = "", gold_struc = "", gold_paired_pos_tuple = None, gold_unpaired_pos = None, pred_paired_pos_tuple = None, pred_unpaired_pos = None, allow_slip=False):
    """
    Evaluates the predicted RNA secondary structure against a gold structure.

    Args:
        pred_struc (str): A string representing the predicted RNA secondary structure, where '(' and ')' denote paired
                          nucleotides and '.' denotes unpaired nucleotides.
        gold_struc (str): A string representing the gold RNA secondary structure, where '(' and ')' denote paired
                          nucleotides and '.' denotes unpaired nucleotides.
        allow_slip (bool): A boolean indicating whether to allow one-nucleotide slips in the predicted structure when
                           calculating precision and sensitivity. If True, a predicted paired nucleotide position is
                           considered correct if it is within one nucleotide of a gold paired nucleotide position.

    Returns:
        tuple: A tuple containing four float values:
            - Precision: The fraction of predicted paired nucleotide positions that are correct.
            - Sensitivity: The fraction of gold paired nucleotide positions that are correctly predicted.
            - F1 score: The harmonic mean of precision and sensitivity.
            - Structural distance: The difference between the length of the gold structure string and twice the number
                                   of common paired nucleotide positions plus the number of common unpaired nucleotide
                                   positions. This represents how close the predicted structure is to the gold structure.

    Raises:
        AssertionError: If the length of the predicted and gold structure strings are not equal.

    Example:
        >>> evaluate('((..((...)).))', '((..((...)).))')
        (1.0, 1.0, 1.0, 0)
        >>> evaluate('((..((...)).))', '((..(.()..).))')
        (0.75, 0.75, 0.75, 4)
    """

    if gold_paired_pos_tuple is None or gold_unpaired_pos is None:
        gold_paired_pos_tuple, gold_unpaired_pos = parse_secondary_structure(gold_struc)
    if pred_paired_pos_tuple is None or pred_unpaired_pos is None:
        pred_paired_pos_tuple, pred_unpaired_pos = parse_secondary_structure(pred_struc)

    
    gold_len = 2 * len(gold_paired_pos_tuple) + len(gold_unpaired_pos) if gold_struc == "" else len(gold_struc)
    assert len(pred_struc) == gold_len, "Length of predicted and gold structure strings must be equal\nPredicted Length: {}\nGold Length: {}".format(
        len(pred_struc), gold_len
    )

    if allow_slip:
        graph = defaultdict(set)
        for (i, j) in pred_paired_pos_tuple:
            for (x,y) in [(i, j), (i-1,j), (i+1,j), (i,j-1), (i,j+1)]:
                if (x, y) in gold_paired_pos_tuple:
                    graph[(i,j)].add((str(x), str(y)))
        
        matching = HopcroftKarp(graph).maximum_matching()
        # only select values that are tuple of string
        common_paired = set([(int(i), int(j)) for (i, j) in matching.values() if isinstance(i, str)])

        all_paired_pos = set()
        pred_new_unpaired_pos = []
        for (i, j) in common_paired:
            all_paired_pos.add(i)
            all_paired_pos.add(j)
        for (i, j) in pred_paired_pos_tuple:
            if (i, j) not in matching:
                all_paired_pos.add(i)
                all_paired_pos.add(j)

        # get new unpaired pos
        for i in range(len(pred_struc)):
            if i not in all_paired_pos and pred_struc[i] != '*':
                pred_new_unpaired_pos.append(i)
        # get common unpaired pos
        common_unpaired = set(pred_new_unpaired_pos).intersection(gold_unpaired_pos)
    else:
        common_paired = set(pred_paired_pos_tuple).intersection(gold_paired_pos_tuple)
        common_unpaired = set(pred_unpaired_pos).intersection(gold_unpaired_pos)


    precision = len(common_paired) / (len(pred_paired_pos_tuple) + 1e-10)
    sensitivity = len(common_paired) / (len(gold_paired_pos_tuple) + 1e-10)
    f1 = 2 * precision * sensitivity / (precision + sensitivity + 1e-10)
    structural_distance = gold_len - (2 * len(common_paired) + len(common_unpaired))
    return precision, sensitivity, f1, structural_distance


def get_alignment_to_sequence_mapping(aligned_sequence):
    """
    Returns a dictionary mapping the indices of the aligned sequence to the indices of the unaligned sequence.

    Args:
        aligned_sequence (str): A string representing the aligned sequence, where '-' denotes a gap.

    Returns:
        dict: A dictionary mapping the indices of the aligned sequence to the indices of the unaligned sequence.

    Example:
        >>> get_alignment_to_sequence_mapping('AUCG-AUCG')
        {0: 0, 1: 1, 2: 2, 3: 3, 4: 3, 5: 4, 6: 5, 7: 6}
    """
    mapping = {}
    j = 0
    for i, x in enumerate(aligned_sequence):
        if x != "-":
            mapping[i] = j
            j += 1
        else:
            mapping[i] = j - 1
    return mapping

def pairs_to_struc(pairs, seq_length):
    struc = ['.' for _ in range(seq_length)]
    for i, j in pairs:
        struc[i] = '('
        struc[j] = ')'
    return ''.join(struc)

def map_consns_struc_to_aln_seq(consns_struc, aln_seq):
    a2s_map = get_alignment_to_sequence_mapping(aln_seq)
    seq = aln_seq.replace("-", "")  # ungapped sequence

    # get the structure corresponding to the unaligned sequence
    struc = ["."] * len(seq)

    for i, j, (b1, b2) in parse_secondary_structure(consns_struc, True)[0]:
        if aln_seq[i] + aln_seq[j] not in VALID_PAIRS:
            continue
        if seq[a2s_map[i]] + seq[a2s_map[j]] not in VALID_PAIRS:
            continue

        struc[a2s_map[i]] = b1
        struc[a2s_map[j]] = b2

    return "".join(struc), seq

def map_consns_bpp_to_aln_seq(consns_bpp, aln_seq, threshold=0.001):
    a2s_map = get_alignment_to_sequence_mapping(aln_seq)
    seq = aln_seq.replace("-", "")  # ungapped sequence

    # get the structure corresponding to the unaligned sequence
    bpp = defaultdict(lambda: defaultdict(float))

    for i in range(len(aln_seq)):
        for j in consns_bpp[i].keys():
            if aln_seq[i] + aln_seq[j] not in VALID_PAIRS:
                continue
            if seq[a2s_map[i]] + seq[a2s_map[j]] not in VALID_PAIRS:
                continue
            if consns_bpp[i][j] >= threshold:
                bpp[a2s_map[i]][a2s_map[j]] = consns_bpp[i][j]
    
    return bpp, seq

def get_struc_from_file(file_path, backward_search=False):
    """
    Get the structure string from a file.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in (lines[::-1] if backward_search else lines):
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] in ['.', '(']:
                return line
            
def get_seq_from_file(file_path, backward_search=False):
    """
    Get the sequence string from a file.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate((lines[::-1] if backward_search else lines)):
            header = lines[i-1].strip() if i > 0 else ""
            line = line.strip()
            if line[0] in set(['A', 'U', 'C', 'G', '-']):
                return header, line
    
def get_bpp_from_file(file_path, threshold=0.01, input_is_zero_indexed=False):
    """
    Get the bpp matrix from a file.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    bpp_matrix = defaultdict(lambda: defaultdict(float))
    seq_length = -1
    for line in lines:
        split = line.split(' ')
        if len(split) == 3 and len(split[0]) < 100 and float(split[2]) >= threshold:
            if input_is_zero_indexed:
                i, j, prob = int(split[0]), int(split[1]), float(split[2])
            else:
                i, j, prob = int(split[0]) - 1, int(split[1]) - 1, float(split[2])
            bpp_matrix[i][j] = prob
            seq_length = max(seq_length, j + 1)
    return bpp_matrix, seq_length

def convert_bpp_to_list(bpp_matrix, seq_length, threshold=0.01):
    consns_bpp = []
    for i in range(seq_length):
        for j in bpp_matrix[i].keys():
            if bpp_matrix[i][j] >= threshold:
                consns_bpp.append((i, j, bpp_matrix[i][j]))
    return consns_bpp


def get_unpaired_probs(paired_probs, seq_length):
    unpaired_probs = {}

    for i in range(seq_length):
        unpaired_probs[i] = 1.00

    for l in paired_probs.keys():
        for r, p in paired_probs[l].items():
            unpaired_probs[l] -= p
            unpaired_probs[r] -= p

    return unpaired_probs


def parse_ct_file(ct_file_path):
    paired_pos = set()
    unpaired_pos = set()
    
    with open(ct_file_path, 'r') as file:
        # Skip header lines if necessary
        next(file) 
        
        for line in file.readlines():
            parts = line.strip().split()
            if len(parts) < 6:
                continue  # skip malformed lines
            
            i = int(parts[0]) - 1   # Adjust index for 0-based Python indexing
            j = int(parts[4]) - 1   # Adjust index for 0-based Python indexing

            if j != -1 and i < j:
                paired_pos.add((i, j))
            elif j == -1:
                unpaired_pos.add(i)

    return paired_pos, unpaired_pos

def calculate_pair_sequence_identity(seq1, seq2):
    """
    Calculate the pairwise sequence identity between two aligned sequences.
    Skips positions where both sequences have gaps.
    
    Args:
    seq1 (str): First aligned sequence with possible gaps ('-')
    seq2 (str): Second aligned sequence with possible gaps ('-')
    
    Returns:
    float: Pairwise sequence identity as a percentage
    """

    # Ensure the sequences are of the same length
    if len(seq1) != len(seq2):
        raise ValueError("The two sequences must be of the same length for comparison, but found {} {} {}".format(seq1, "\n\n", seq2))
    
    # Initialize counters
    identical_positions = 0
    aligned_positions = 0
    
    # Loop through each position to compare the two sequences
    for base1, base2 in zip(seq1, seq2):
        # Skip positions where both sequences have gaps
        if base1 == '-' and base2 == '-':
            continue
        
        # Increment the counter for aligned positions
        aligned_positions += 1
        
        # Count identical positions
        if base1 == base2:
            identical_positions += 1
    
    # If no aligned positions are present, return 0 to avoid division by zero
    if aligned_positions == 0:
        return 0.0
    
    # Calculate the pairwise sequence identity as a percentage
    return (identical_positions / aligned_positions) * 100

def calculate_msa_seq_identity(msa):
    """
    Calculate the average pairwise sequence identity for a multiple sequence alignment (MSA).
    
    Args:
    msa (list of str): List containing aligned sequences.
    
    Returns:
    float: Average pairwise sequence identity.
    """
    num_sequences = len(msa)
    if num_sequences < 2:
        raise ValueError("At least two sequences are required to compute pairwise sequence identity")
    
    total_identity = 0
    num_comparisons = 0
    
    # Calculate pairwise identity for each unique pair of sequences
    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            identity = calculate_pair_sequence_identity(msa[i], msa[j])
            total_identity += identity
            num_comparisons += 1
    
    # Calculate the average identity
    average_identity = total_identity / num_comparisons if num_comparisons > 0 else 0
    return average_identity
