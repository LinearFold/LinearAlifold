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

- `-b`, `--b` `<int>`: Set the beam size. Default is 100.
- `-pt`, `--partition`: Enable the partition mode.
- `--verbose`: Enable verbose mode to print out more information. By default, this is set to false.
- `--em` `<int>`: Select the energy model. Choose between 1 (Vienna) or 2 (BL*). Default is 2.
- `-ct`, `--cutoff` `<int>`: Set the conservation score cutoff threshold. Default is -40.
- `-bt`, `--beta` `<float>`: Specify the beta value. Default is 1.2.
- `-dt`, `--delta` `<float>`: Specify the delta value. Default is 0.1.
- `-tt`, `--threshknot-threshold` `<float>`: Set the threshknot threshold. Default is 0.3.
- `--bpp-file` `<string>`: Path to the bpp file. When it is specified, the program will output the base pair probability matrix to the given file.
- `--mea-file` `<string>`: Path to the mea file. When it is specified, the program will output the MEA structure to the given file.
- `--threshknot-file` `<string>`: Path to the threshknot file. When it is specified, the program will output the ThreshKnot structure to the given file.

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
