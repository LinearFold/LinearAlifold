## Issues with the Old Codebase 
The entire codebase is excessively verbose and poorly structured. The excessive redundancy and logic misplacements within the codebase, lead to an error prone and challenging-to-maintain code environment.

### Redundant Variables Leading to Clutter and Confusion
- The codebase includes numerous variables which, although used, contribute to redundancy. This creates a cluttered and confusing layout, making the code difficult to read and maintain. Notable examples include:

  - **Redundant Variables in Loops:** There is an `unordered_map<int, pair<int, State>>` named `bestP`, yet it redundantly uses `beamstepP = bestP[j]` within a for-loop, adding no extra functionality. Similar redundancy is observed with `bestM`, `bestMulti`, `bestM2`, `bestH`, and `bestC`.

  - **Redundant Index Variables:** The `a2s_fast` (alignment to sequence index) variable is present, but redundant variables such as `a2s_i`, `a2s_j` (equivalent to `a2s_fast[i]` and `a2s_fast[j]`) are declared multiple times in different beams. This pattern is repeated for `a2s_i_1` (=`a2s_fast[i-1]`), `a2s_fast_j_1` (=`a2s_fast[j-1]`), and other similar variables like `SS_fast_i`, `SS_fast_i_1`, `s3_fast_i`, `s5_fast_j`.

  - **Redundant Data Structure Initialization:** The `next_pair` data structure is initialized redundantly and serves no practical use, as its functionality is already covered by `next_pair_msa`. This unnecessary initialization of `next_pair` adds to the clutter and serves no purpose.

### Lack of Modularization

- **Repeated Logic Instead of Reusable Functions:**
  - The codebase demonstrates a significant lack of modularization, notably evident in the repetitive implementation of logic for various functionalities like **jump next**. This leads to an unnecessary increase in the volume of code and complicates maintenance and updates.

  - **Specific Example - Repetition in Calculating `next_position[i][j]`:** A specific instance of this issue is the repeated logic for calculating the value of `next_position[i][j]` across different parts of the code. Instead of centralizing this calculation in a single, reusable function, the codebase redundantly implements the same logic in multiple locations. This not only bloats the code but also poses a high risk for errors and inconsistencies, especially when modifications or bug fixes are needed.

- **Improper Utilization of Initialization Functions:**
  - Despite having a designated `prepare` function for initializing all necessary data structures prior to running the parsing function, the codebase inconsistently handles some initializations.

  - **Initialization of `next_position` and `next_pair_MSA`:** Notably, the initialization of `next_position` and `next_pair_MSA` occurs within the parsing function, despite the existence of a separate `prepare` function meant for such initializations. Ideally, all initial setups, including those for `next_position` and `next_pair_MSA`, should occur within the `prepare` function to maintain code cleanliness and clarity.

### Misplaced Initializations in Separate Functionality File
- **Inconsistent Initialization Practices:**
  - Although there is a designated `prepare` function for all initializations, certain data structures are incorrectly initialized in a separate functionality file, `ribo.h`, which serves a completely different purpose. This file is responsible for the initialization of:
    - `vector<vector<int>> a2s_fast`
    - `vector<vector<int>> s5_fast`
    - `vector<vector<int>> s3_fast`
    - `vector<vector<int>> SS_fast`
    - `vector<float> smart_gap`
  
  - This separation of initializations from the `prepare` function disrupts the logical structure of the code and contributes to confusion. 

## Refactoring Summary
These refactoring efforts significantly improved the structure, clarity, and maintainability of the codebase, resolving existing problems and laying a solid foundation for future development.

### Unified MFE and Partition Algorithm into a Single Program
- Previously, the MFE (Minimum Free Energy) and partition algorithm were implemented as separate programs, each with its own header file. 
- Simplified this by creating a single program with two distinct files: `inside.cpp` and `outside.cpp`. This change not only reduced redundancy but also brought about greater coherence in the codebase. The MFE and the inside algorithm, being very similar except for their summary operators, now share the same `inside.cpp` file, demonstrating an efficient use of shared logic and resources.

### Rewrite and Improved Structuring
- Removed redundant variables and inconsistencies, thereby optimizing the code for improved performance.
- Rewrote the entire inside and outside parts from scratch to ensure a modularized, well-structured, and a more readable codebase. This guaranteed error-free logic in the code and enhanced functionality and performance.
- Organized the codebase into separate files, achieving a cleaner and more logical structure. Created separate files like `linearalifold.h`, `linearalifold.cpp`, `inside.cpp`, `outside.cpp`.

### Proper Use of Functions and Modularization
- Correctly utilized functions, such as employing the `prepare` function for all initializations, which enhanced the code's organization.
- Implemented effective modularization by creating reusable functions for repeated logic, simplifying the code and making it easier to manage and debug.

## Changes Made and New Features
C++ is a sensitive language (especially when using pointers) where errors can easily arise from verbose, redundant, and poorly structured code. Due to the comprehensive rewrite of the major parts of the codebase from scratch, it was challenging to identify every inconsistency or error in the original code. The following sections detail the major changes, inconsistencies, and errors identified and addressed during the refactoring process.

- **Beam of H:**
  - Implemented a correct constraint for Hairpins in the new codebase with the following logic:
    ```cpp
    int jnext = get_cache_next_position(j, j + 3);
    while (jnext != -1 && (smart_gap[jnext] - smart_gap[j] < 4 * smart_gap_threshold)) {
        jnext = get_cache_next_position(j, jnext);
    }
    ```
  - Removed unnecessary `if` condition from the old codebase as seen in this [code segment](https://github.com/LinearFold/LinearAlifold/blob/d801ac370677ed24ec9634567cf9a234bd92c630/LinearAlifold_MFE/src/Linearalifold.cpp#L408).

- **Beam of M2:**
  - Revised the max length condition for M2 -> Multi jump:
    - Old version: 
      ```cpp
      for (int p = i-1; ((smart_i - smart_gap[p]) <= 2*SINGLE_MAX_LEN) and p >= 0; --p)
      ```
    - New version: 
      ```cpp
      for (int p = i - 1; p >= 0 && (smart_gap[i - 1] - smart_gap[p] <= SINGLE_MAX_LEN); --p)
      ``` 
  - Removed unnecessary `if` condition from the old codebase as seen in this [code segment](https://github.com/LinearFold/LinearAlifold/blob/d801ac370677ed24ec9634567cf9a234bd92c630/LinearAlifold_MFE/src/Linearalifold.cpp#L1014).

- **Cube Pruning:**
  - If compiled with the cube pruning flag, the old version only sorted `bestM` in the beam of M1. It lacked implementation for M2 = M + P with cube pruning.
  - The new version correctly employs cube pruning (inherited from LinearFold).

- **Correct Jump Next:**
  - The old codebase lacked a separate function for the jump next feature, leading to repeated computations across multiple beams. This resulted in discrepancies between the old and new versions, with the new version being more simple and correct.
  - The new codebase introduces separate functions for jump next and jump previous (for lazy outside), enhancing readability, modularization, and maintainability.

- **C = C + P Rule:**
  - The old codebase implemented an inconsistent check for identifying the first non-gap nucleotide to the left of index `i`, and to the right of index `j`.

  - **Old Code:**
    - In the old version, the logic used can be seen in this [code segment](https://github.com/LinearFold/LinearAlifold/blob/d801ac370677ed24ec9634567cf9a234bd92c630/LinearAlifold_MFE/src/Linearalifold.cpp#L818):
      ```cpp
      int new_nuck = (a2s_i[s] > 0) ? s5_i[s] : -1;
      int new_nucj1 = (a2s_j[s] < a2s_seq_length_1[s]) ? s3_j[s] : -1;
      ```
  - **New Code:**
    - The new version employs a more consistent approach:
      ```cpp
      int nuci_1 = (i - 1) > -1 ? s5_fast[i][s] : -1;
      int nucj1 = (j + 1) < seq_length ? s3_fast[j][s] : -1;
      ```

- **Beam of P:**
  - Implemented a revised condition for new helix/single branch loop in the Beam of P using `smart_gap`.

  - **New Version:**
    ```cpp
    for (int p = i - 1; p >= 0 && (smart_gap[i - 1] - smart_gap[p] <= SINGLE_MAX_LEN); --p) {
        int q = get_cache_next_position(p, j);
        while (q != -1 && ((smart_gap[i - 1] - smart_gap[p]) + (smart_gap[q - 1] - smart_gap[j]) <= SINGLE_MAX_LEN)) {
            // new state is of shape p..i..j..q
            ...
            ...
        }
    }
    ```

  - **Old Version:**
    ```cpp
    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
    // block of jump next code
    ...
    ...
        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
            // new state is of shape p..i..j..q
            ...
            ...
        }
    }
    ```

- **Incorrect Pair Score Computation Correction:**
  - The old version had an issue with the incorrect computation of the pair score, specifically in the starting value of the inner loop variable `l` and not averaging the final score by dividing by the number of sequences.

  - **Old Version: (Wrong)**
    ```cpp
    double score = 0.;
    for (int k = 1; k <= 6; k++) { /* ignore pairtype 7 (gap-gap) */
        for (int l = k; l <= 6; l++) {
            score += (double)pfreq[k] * (double)pfreq[l] * (double)ribo[k][l];
        }
    }
    return 100 * score / n_seq - 100 * (pfreq[0] + pfreq[7] * 0.25);
    ```

  - **New Version: (Correct)**
    ```cpp
    double score = 0.;
    for (int k = 1; k <= 6; k++) { /* ignore pairtype 7 (gap-gap) */
        score += (double)pfreq[k] * (double)(pfreq[k] - 1) / 2 * (double)ribo[k][k];
        for (int l = k + 1; l <= 6; l++) {
            score += (double)pfreq[k] * (double)pfreq[l] * (double)ribo[k][l];
        }
    }
    return beta * (100 * score / n_seq - delta * 100 * (pfreq[0] + pfreq[7] * 0.25)) / n_seq;
    ```
  - Additionally, the new version incorporates parameters `beta` and `delta` for further tuning of the covariance term and the penalty term.


- **MFE Score Averaging Correction:**
  - The old version of the code did not average the energies of different loops (by dividing by number of sequences), which was an incorrect approach. The new version corrects this issue by averaging the energies, ensuring a more accurate and representative Minimum Free Energy (MFE) score calculation.

- **Optimization in Division Operations:**
  - The new version optimizes computational efficiency by replacing division operations with multiplication. This is achieved by pre-calculating the inverse of certain variables and using these inverses for multiplication, which is computationally faster. 
    - In the old version, division operations like `score / n_seq` or `score / kTn` were used.
    - In the new version, these are replaced by:
      - Pre-calculating `inv_n = 1 / n_seq` and `inv_ktn = 1 / kTN`.
      - Then applying multiplication: `score *= inv_n` or `score *= inv_ktn`.

- **Unpaired bases in  Multi Loop:**
  - The old version neglected to count the number of unpaired bases in Multi Loops due to the energy being set to 0 in the Vienna energy model. In contrast, the new version supports a detailed calculation for these unpaired bases. It also includes an approximate calculation method similar to that used in LinAliFold.

- **Implementation of Multi Jump Rule with Smart Gap:**
  - In the new version, the multi jump rule is implemented using `smart_gap` and is restricted to a maximum jump length of 30.
  - **Code Snippet:**
    ```cpp
    // 2. extend (i, j) to (i, jnext)
    int jnext = get_cache_next_position(i, j);
    if (jnext != -1 && (smart_gap[jnext - 1] - smart_gap[j - 1] <= SINGLE_MAX_LEN)) {
        ...
        ...
    }
    ```

- **Support for Multiple Energy Models:**
  - Whereas the old version was limited to supporting only the Vienna energy model, the new version expands this capability. It now supports both the Vienna and BL* energy models. Additionally, it introduces the functionality to parse a new energy model text file for parameters.

- **Faster Lazy Outside Algorithm:**
  - The new version introduces a significantly faster variant of the outside algorithm, called 'Lazy Outside.' This updated algorithm completes in approximately 1-2% of the time taken by the original outside algorithm.

- [TODO] **Centroid Structure Computation:**
    - The new version will include support for computing Centroid Structures in addition to MFE and MEA structures.

- [TODO] **More Efficient Backtracking Algorithm:**
    - The new version will enhance the backtracking algorithm, aiming to make it more efficient. This improvement will further optimize the overall performance of the software.
