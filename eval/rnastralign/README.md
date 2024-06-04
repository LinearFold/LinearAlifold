# RNAstralign Evaluation
This directory contains the structural distance evaluation script and RNAstralign dataset.

## Data
Four families (Group I Intron, tmRNA, tRNA, and 5S rRNA) are used for parameter tuning and another four families (SRP, RNaseP, telomerase, and 16S rRNA) are used for testing. For Group I Intron, 5S rRNA, SRP, RNaseP, and 16S rRNA, there are multiple subfamilies within each family, so we chose one specific subfamily for these five families (See table below for more details).

| family     | subfamily           | avg. seq. len. | avg. seq. identity |
|------------|---------------------|----------------|--------------------|
| Group 1    | IC1                 |          428.5 |               0.31 |
| tmRNA      | -                   |          367.4 |               0.35 |
| tRNA       | -                   |          77.1  |               0.48 |
| 5S rRNA    | Bacteria            |          116.2 |               0.61 |
|------------|---------------------|----------------|--------------------|
| SRP        | Protozoan           |          285.8 |               0.35 |
| RNaseP     | Bacterial           |          360.0 |               0.43 |
| telomerase | -                   |          444.9 |               0.45 |
| 16S RNA    | Alphaproteobacteria |         1419.2 |               0.85 |

There are two versions of the data, aligned version (all the homologs in the sample are aligned) and unaligned version:
- Aligned Version: [data/aln/](./data/aln/)
- Unaligned Version: [data/no_aln/](./data/no_aln/)

