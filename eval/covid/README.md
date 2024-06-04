# Covid Evaluation
This directory contains the evaluation scripts (Structural Distance and Ensemble Defect) and data for SARS-CoV-2 and diverse SARS-related genomes.

## Data
The dataset contains the reference sequence, SARS-CoV-2 variants (Alpha, Beta, Delta, and Omicron), and SARS-related genomes.
The table below shows the number of genomes for each $k$ (number of sequences in the MSA). You can find the dataset in the [/data/v1/](./covid/data/v1/) directory.

**Table: SARS-CoV-2 and SARS-related datasets. Ref is the SARS-CoV-2 reference sequence, Alpha–Delta are the SARS-CoV-2 variants, and SARSr are SARS-related genomes.**
| k   | Ref | Alpha | Beta | Delta | Omicron | SARSr |
|-----|-----|-------|------|-------|---------|-------|
| 10  | 1   | 2     | 2    | 1     | 1       | 3     |
| 30  | 1   | 4     | 5    | 4     | 4       | 12    |
| 50  | 1   | 7     | 8    | 7     | 7       | 20    |
| 100 | 1   | 14    | 15   | 15    | 15      | 40    |
| 200 | 1   | 33    | 33   | 33    | 20      | 80    |
| 300 | 1   | 59    | 60   | 40    | 20      | 120   |
| 400 | 1   | 79    | 80   | 60    | 20      | 160   |

## Hybrid Reference Structure
To get the hybrid reference structure, we combined the experimentally guided structures from Huston et al. and the experimentally determined end-to-end pairs (Arch3, ranges from (60,29868) to (80,29847)) from Ziv et al. by the following steps:

1. Get (local) structures in 5’ and 3’ UTR regions from Huston et al. (the 5’ UTR ranges from 1 to 400 and the 3’ UTR from 29543 to 29876 on the reference sequence).
2. Remove (local) pairs (i,j) from the structures if i or j is in the global Arch3 pairs. These local pairs were predicted by the local folding software which can only predict pairs within a local window.
3. Combine the modified structures and the end-to-end Arch3 pairs from Ziv et al.

The dashes ('$-$') in the sequence and structure below represent the split point between the 5' and 3' regions.
```
>NC_045512.2_Wuhan_seafood_market_pneumonia_virus_isolate_Wuhan-Hu-1__complete_genome
AUUAAAGGUUUAUACCUUCCCAGGUAACAAACCAACCAACUUUCGAUCUCUUGUAGAUCUGUUCUCUAAACGAACUUUAAAAUCUGUGUGGCUGUCACUCGGCUGCAUGCUUAGUGCACUCACGCAGUAUAAUUAAUAACUAAUUACUGUCGUUGACAGGACACGAGUAACUCGUCUAUCUUCUGCAGGCUGCUUACGGUUUCGUCCGUGUUGCAGCCGAUCAUCAGCACAUCUAGGUUUCGUCCGGGUGUGACCGAAAGGUAAGAUGGAGAGCCUUGUCCCUGGUUUCAACGAGAAAACACACGUCCAACUCAGUUUGCCUGUUUUACAGGUUCGCGACGUGCUCGUACGUGGCUUUGGAGACUCCGUGGAGGAGGUCUUAUCAGAGGCACGUCAACAU----------ACUCAACUCAGGCCUAAACUCAUGCAGACCACACAAGGCAGAUGGGCUAUAUAAACGUUUUCGCUUUUCCGUUUACGAUAUAUAGUCUACUCUUGUGCAGAAUGAAUUCUCGUAACUACAUAGCACAAGUAGAUGUAGUUAACUUUAAUCUCACAUAGCAAUCUUUAAUCAGUGUGUAACAUUAGGGAGGACUUGAAAGAGCCACCACAUUUUCACCGAGGCCACGCGGAGUACGAUCGAGUGUACAGUGAACAAUGCUAGGGAGAGCUGCCUAUAUGGAAGAGCCCUAAUGUGUAAAAUUAAUUUUAGUAGUGCUAUCCCCAUGUGAUUUUAAUAGCUUCUUAGGAGAAUGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
......(((((.(((((....)))))..)))))...........(((((.....)))))((((((((.((.(((((((((...((((((((.((.((((.(((.....))).)))))).))))))))......((((.....))))...(((((((((((..(((((...(((.(((((((((((..((((((.(((((......)))))..))))))......)))(((((((.((......)))))))))(((....)))))))))))))).))))).))))...))))))).......((((((...........((((((...))))))....)))))).....(((((.(((((((((((((.....)))).))))..))))).)))))......----------..........((((((.((..((((.......))))..))...)))))).......((((((..((.(((((((((((..(((...((......))...))).)))))))))))))))))))..............((((((((((...........))))))))).).......((((((((((......(((.........(((((((((............(((..((.(((((((((....((.((...((......))..))))...))))).))))...))....)))...............))))))))).........))).......))))))......))))..)))).)).))))).))))))))..........
```

Here's the hybrid reference structure text file (dot-bracket): [hybrid_reference_structure.fasta](./data/hybrid_reference_structure.fasta).\
We have also included the drawing/visualization of the hybrid reference structure in .nsd format here: [hybrid_reference_structure.nsd](./data/hybrid_reference_structure.nsd). This file can be opened using the RNAstructure Editor software (link: https://rna.urmc.rochester.edu/RNAstructureDownload.html).