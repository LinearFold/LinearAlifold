# LinAliFold and CentroidLinAliFold
LinAliFold and CentroidLinAliFold is fast MFE and MEA-based consensus secondary structure prediction software from RNA alignments using beam search, respectively.

## Usage
Both tools can take the multi-FASTA format as an input file format, and predict its common secondary structure.

An example:
```
$ ./bin/LinAliFold -i sample.fasta
score:-23.13
(((((((..((((...........)))).(((((.......)))))...................(((((.......))))))))))))....
```
## Options   
    LinAliFold Options:
    (Required)
        -i STR    InputFileName
        
    (Optional) 
        -b INT    The beam size [default:100]
        -d INT    The value of the parameter delta [default: 1.0]
        -e INT    The value of the parameter beta [default: 1.0]
        -r INT    Designation of the covariation score model 0: simple score model, 1: ribosum score model [defualt:0]
                        
    CentroidLinAliFold Options:
    (Required)
        -i STR    InputFileName
        
    (Optional)
        -b INT    The beam size [default:100]
        -d INT    The value of the parameter delta [default: 1.0]
        -e INT    The value of the parameter beta [default: 1.0]
        -r INT    Designation of the covariation score model 0: simple score model, 1: ribosum score model [defualt:0]
        -g INT    The value of the parameter gamma [default: 1.0]
        -t DBL    The threshold value for base pairing probabilities [default: 0.0001]
        -w DBL    The value of the mixture weight parameter [default: 0.5]
        -p INT    Designation of the selection method of gamma 0: user-defined parameter 1: maximization of the pseudo-expected accuracy of MCC [default: 0]
        -o INT    Designation of the output style 0: the bpp is not outputted 1: the bpp is outputted [default: 0]

For the meaning of parameters delta and beta, please see the following RNAaliFold paper [1].

[1] Bernhart, S.H., Hofacker, I.L., Will, S. et al. RNAalifold: improved consensus structure prediction for RNA alignments. BMC Bioinformatics 9, 474 (2008)

## Version
Version 1.0.0 (2022/06/17)

## Acknowledgements
In developing this program, we referred to the source code from the LinearFold and LinearPartition. We thank Dr. Liang Huang and LinearFold development group. You can download the source code of LinearFold and LinearPartition from https://github.com/LinearFold/LinearFold and https://github.com/LinearFold/LinearPartition. 

## License
This software is released under the MIT License, see LICENSE.txt.  

## Reference
Tsukasa Fukunaga and Michiaki Hamada. "LinAliFold and CentroidLinAliFold: Fast RNA consensus secondary structure prediction for aligned sequences using beam search methods." under submission.
