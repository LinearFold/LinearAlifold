def parse_secondary_structure(struc):
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
                paired_indices.append((stack1.pop(), i))
            elif x == "]":
                paired_indices.append((stack2.pop(), i))
            elif x == "}":
                paired_indices.append((stack3.pop(), i))
            elif x == ">":
                paired_indices.append((stack4.pop(), i))
            elif x == ".":
                unpaired_indices.append(i)
    except Exception as _:
        print("[Error] Unbalanced parenthesis in structure string")
        return None

    if stack1 or stack2 or stack3 or stack4:
        print("[Error] Unbalanced parenthesis in structure string")

    return paired_indices, unpaired_indices


def evaluate(pred_struc, gold_struc, allow_slip=False):
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
    assert len(pred_struc) == len(
        gold_struc
    ), "Length of predicted and gold structure strings must be equal\nPredicted Length: {}\nGold Length: {}".format(
        len(pred_struc), len(gold_struc)
    )

    pred_paired_pos, pred_unpaired_pos = parse_secondary_structure(pred_struc)
    gold_paired_pos, gold_unpaired_pos = parse_secondary_structure(gold_struc)

    common_unpaired = set(pred_unpaired_pos).intersection(gold_unpaired_pos)
    if allow_slip:
        common_paired = []
        for i, j in pred_paired_pos:
            if set([(i, j), (i, j - 1), (i, j + 1), (i - 1, j), (i + 1, j)]) & set(gold_paired_pos):
                common_paired.append((i, j))
        common_paired = set(common_paired)
    else:
        common_paired = set(pred_paired_pos).intersection(gold_paired_pos)

    valid = True


    precision = len(common_paired) / (len(pred_paired_pos) + 0.000001)
    sensitivity = len(common_paired) / (len(gold_paired_pos) + 0.000001)
    f1 = 0 if precision * sensitivity == 0 else 2 * precision * sensitivity / (precision + sensitivity)
    structural_distance = len(gold_struc) - (2 * len(common_paired) + len(common_unpaired))

    if (f1 > 1):
        f1 = 1
    

    return valid, precision, sensitivity, f1, structural_distance


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
