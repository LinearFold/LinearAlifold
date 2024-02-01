#!/bin/bash
# cat ../eval/Data/samples/sample07.fasta | ./bin/laf_mfe_vienna 100 1 1 1 -40 1 1 1


cat ./eval/samples/sample17.fasta | valgrind --leak-check=full --track-origins=yes -s ./bin/laf_mfe_vienna 100 1 1 1 -40 1 1 0 "" ""