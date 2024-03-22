#!/usr/bin/env python3

import os
import subprocess
import sys
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="LinearAlifold: Linear-Time Consensus Structure Prediction for RNA Alignments")

    parser.add_argument("-v", "--verbose", action="store_true", default=False, help="Enable verbose mode to print out detailed information (DEFAULT=FALSE)")
    parser.add_argument("-b", "--beam", type=int, default=100, help="Set beam size (beam=0 will result in exact search) (DEFAULT=100)")
    parser.add_argument(
        "-e", "--energy-model", type=int, choices=[1, 2], default=2, help="Energy model choice: 1 (Vienna) or 2 (BL*) (DEFAULT=2)"
    )
    parser.add_argument('-p', '--partition', action='store_true', help="Partition Mode: computes BPPs, MEA Structure, and Threshknot structure")
    parser.add_argument('-s', '--sample', type=int, default=0, help="Sampling Mode: sample size for running sampling mode (0 means sampling mode is off) (DEFAULT=0)")

    parser.add_argument("-c", "--cutoff", type=int, default=-40, help="Set pairability/conservation score cutoff threshold (DEFAULT=-40)")
    parser.add_argument("-y", "--beta", type = float, default=1.2, help="Set parameter beta value (DEFAULT=1.2)")
    parser.add_argument("-z", "--delta", type = float, default=0.1, help="Set parameter delta value (DEFAULT=0.1)")
    parser.add_argument("-g", "--gamma", type=float, default=2, help="Set parameter gamma value (DEFAULT=0)")
    parser.add_argument('-t', '--threshknot-threshold', type=float, default=0.3, help='Set threshknot threshold')
    
    parser.add_argument('-B', "--bpp-file", type=str, help="BPP matrix save path", default="")
    parser.add_argument('-M', "--mea-file", type=str, help="MEA structure save path", default="")
    parser.add_argument('-T', "--threshknot-file", type=str, help="Threshknot structure save path", default="")
    parser.add_argument('-C', "--centroid-file", type=str, help="Centroid structure save path", default="")


    parser.add_argument('-m', "--multi-approx", action="store_true", help="Use multiloop approximation (DEFAULT=FALSE)")
    parser.add_argument('-o', "--non-lazy", action="store_false", help="Do not use lazy outside (DEFAULT=TRUE)")
    return parser.parse_args()

def main():
    args = get_args()

    options = [str(args.beam),
            "1" if args.verbose else "0",
            str(args.energy_model),
            "1" if args.multi_approx else "0",
            "1" if args.partition else "0",

            str(args.cutoff),
            str(args.beta),
            str(args.delta),
            str(args.gamma),
            
            "1" if args.non_lazy else "0",
            str(args.threshknot_threshold),
            str(args.sample),

            args.bpp_file,
            args.mea_file,
            args.threshknot_file,
            args.centroid_file
        ]
    
    cmd = ["%s/%s" % (os.path.dirname(os.path.abspath(__file__)), \
                      "bin/laf_mfe_vienna" if args.energy_model == 1 else "bin/laf_mfe_bl")] + options
    
    # start the subprocess
    process = subprocess.Popen(cmd, stdin=sys.stdin, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # capture and print the output after the process completes
    stdout, stderr = process.communicate()
    if stdout:
        print(stdout.decode())
    if stderr:
        print("[ERROR] See below for error message")
        print(stderr.decode())

if __name__ == "__main__":
    main()
