#!/usr/bin/env python3

import os
import subprocess
import sys
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="Command-line arguments parser")

    parser.add_argument("-b", "--b", type=int, default=100, help="set beam size, (DEFAULT=100)")
    parser.add_argument("--verbose", action="store_true", default=False, help="print out runtime in seconds, (DEFAULT=FALSE)")

    parser.add_argument(
        "--em", type=int, choices=[1, 2], default=2, help="energy model choice, 1 (Vienna) or 2 (BL*), (DEFAULT=2)"
    )
    parser.add_argument("--multi-approx", type=int, default=0, help="use multiloop approximation, (DEFAULT=0 (False))")
    parser.add_argument("-ct", "--cutoff", type=int, default=-40, help="cutoff threshold (DEFAULT=-40)")
    parser.add_argument("-bt", "--beta", type = float, default=1.2, help="beta value (DEFAULT=1.2)")
    parser.add_argument("-dt", "--delta", type = float, default=0.1, help="delta value (DEFAULT=0.1)")

    parser.add_argument('--bpp-file', type=str, help='bpp file', default="")
    parser.add_argument('--mea-file', type=str, help='mea file', default="")
    parser.add_argument('--threshknot-file', type=str, help='threshknot file', default="")
    parser.add_argument('-tt', '--threshknot-threshold', type=float, default=0.3, help='threshknot threshold')

    parser.add_argument('-pt', '--partition', action='store_true', help='run partition function')
    parser.add_argument('-lz', '--lazy', action='store_true', help='run lazy outside', default=True)
    return parser.parse_args()

def main():
    args = get_args()
    beamsize = str(args.b)
    is_verbose = "1" if args.verbose else "0"
    path = os.path.dirname(os.path.abspath(__file__))
    # print(args)

    options = [beamsize,
            is_verbose,
            str(args.em),
            str(args.multi_approx),
            str(args.cutoff),
            str(args.beta),
            str(args.delta),
            "1" if args.partition else "0",
            "1" if args.lazy else "0",
            str(args.threshknot_threshold),
            args.bpp_file,
            args.mea_file,
            args.threshknot_file,
        ]

    if args.em == 1:
        cmd = [
            "%s/%s" % (os.path.dirname(os.path.abspath(__file__)), "bin/laf_mfe_vienna"),
        ] + options

    elif args.em == 2:
        cmd = [
            "%s/%s" % (os.path.dirname(os.path.abspath(__file__)), "bin/laf_mfe_bl"),
        ] + options

    # Start the subprocess
    process = subprocess.Popen(cmd, stdin=sys.stdin, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Capture and print the output after the process completes
    stdout, stderr = process.communicate()
    if stdout:
        print(stdout.decode())
    if stderr:
        print("There was an error!")
        print(stderr.decode())

if __name__ == "__main__":
    main()
