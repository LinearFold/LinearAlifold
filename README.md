# LinearAlifold: Linear-Time Consensus Structure Prediction for RNA Alignments

Apoorv Malik†, Liang Zhang†, Ning Dai, Sizhen Li, He Zhang, David H. Mathews, Liang Huang†*

\* corresponding author
\† equal contribution

## Description

LinearAlifold is an efficient algorithm for folding aligned RNA homologs that scales linearly with both the sequence length and the number of sequences, based on our work LinearFold that folds a single RNA in linear time. LinearAlifold supports four modes: minimum free energy (MFE), maximum expected accuracy (MEA), ThreshKnot, and stochastic sampling, each of which takes under an hour for hundreds of SARS-CoV variants.

Our work is orders of magnitude faster than RNAalifold (e.g., 0.7 hours on the above 400 genomes, or ∼36× speedup) and achieves higher accuracies when compared to a database of known structures. More interestingly, LinearAlifold’s prediction on SARS-CoV-2 correlates well with experimentally determined structures, outperforming RNAalifold.

## Dependencies
GCC 4.8.5 or above; 
python3

## To Compile
```
make
```

## To Run
(input: a Multiple Sequence Alignment (MSA)):
```
cat MSA_file | ./linearalifold [OPTIONS]
```

## Options   
    LinearAlifold Options:
    (Optional)
        -b,  --beam-size INT                 The beam size [default: 100]
        -vb, --verbose                       Enable verbose mode to print out more information [default: False].
        -em, --energy-model INT              Select the energy model. Choose between 1 (Vienna) or 2 (BL*) [default: 2].
        -ct, --cutoff INT                    Set the conservation score cutoff threshold [default: -40].
        -bt, --beta FLOAT                    Specify the beta value [default: 1.2].
        -dt, --delta FLOAT                   Specify the delta value [default: 0.1].
        -tt, --threshknot-threshold FLOAT    Set the threshknot threshold [default: 0.3].
        -bf, --bpp-file STRING               Path to the bpp file. Outputs the base pair probability matrix if specified.
        -mf, --mea-file STRING               Path to the mea file. Outputs the MEA structure if specified.
        -tf, --threshknot-file STRING        Path to the threshknot file. Outputs the ThreshKnot structure if specified.
        -pt, --partition                     Enable the partition mode.

## Examples

### MFE Mode: Computes the Minimum Free Energy (MFE) Structure
To run LinearAlifold with the MFE mode, you can use the following command:
```
cat ./test_samples/sample01.fasta | ./linearalifold.py 
Minimum Free Energy: -6.82 kcal/mol

MFE Structure: 
...((.(((((((((...........))))))))).)) (-6.82 = -0.72 + -6.09)
```

### Partition Mode: Computes both ThreshKnot and Maximum Expected Accuracy (MEA) Structures
To run LinearAlifold with the ThreshKnot and MEA modes, you can use the following command:
```
cat ./test_samples/sample01.fasta | ./linearalifold.py -pt
Free Energy of Ensemble: -0.06 kcal/mol

MEA Structure:
...((.(((((((((...........))))))))).))


Threshknot Structure:
...((.(((((((((...........))))))))).))
```
