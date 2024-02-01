#!/usr/bin/env python3

import os
import subprocess
import sys
import argparse
import psutil
import time

def get_args():
    parser = argparse.ArgumentParser(description="Command-line arguments parser")

    parser.add_argument("-b", "--b", type=int, default=100, help="set beam size, (DEFAULT=100)")
    parser.add_argument("--verbose", action="store_true", default=True, help="print out runtime in seconds, (DEFAULT=FALSE)")

    parser.add_argument(
        "--em", type=int, choices=[1, 2], default=1, help="energy model choice, 1 (Vienna) or 2 (BL*), (DEFAULT=1)"
    )
    parser.add_argument("--multi-approx", type=int, default=1, help="use multiloop approximation, (DEFAULT=1 (TRUE))")
    parser.add_argument("-ct", "--cutoff", type=int, default=-40, help="cutoff threshold (DEFAULT=-40)")
    parser.add_argument("-bt", "--beta", type = float, default=1.0, help="beta value (DEFAULT=1.0)")
    parser.add_argument("-dt", "--delta", type = float, default=1.0, help="delta value (DEFAULT=1.0)")


    parser.add_argument('--bpp-file', type=str, help='bpp file', default="")
    parser.add_argument('--mea-file', type=str, help='mea file', default="")
    parser.add_argument('--threshknot-file', type=str, help='threshknot file', default="")
    parser.add_argument('-tt', '--threshknot-threshold', type=float, default=0.3, help='threshknot threshold')

    parser.add_argument('-pt', '--partition', action='store_true', help='run partition function')
    parser.add_argument('-lz', '--lazy', action='store_true', help='run lazy outside', default=False)
    return parser.parse_args()

def monitor_memory_usage(process, interval=1):
    start_time = time.time()  # Record the start time
    total_measurements = 0  # Count of measurements
    rolling_average = 0.0  # Rolling average of memory usage
    min_memory = float('inf')  # Initialize min memory to infinity
    max_memory = float('-inf')  # Initialize max memory to negative infinity

    try:
        cpp_process = psutil.Process(process.pid)

        while True:
            if process.poll() is not None:
                break  # Exit the loop if the process has terminated

            new_memory_usage = cpp_process.memory_info().rss / (1024 * 1024)  # Memory usage in MB
            total_measurements += 1
            rolling_average += (new_memory_usage - rolling_average) / total_measurements
            min_memory = min(min_memory, new_memory_usage)
            max_memory = max(max_memory, new_memory_usage)
            print(f"Current Memory Usage: {new_memory_usage:.3f} MB, Rolling Average: {rolling_average:.3f} MB")
            time.sleep(interval)

    except psutil.NoSuchProcess:
        print("Process terminated")
    except Exception as e:
        print(f"An error occurred: {e}")

    end_time = time.time()  # Record the end time
    execution_time = end_time - start_time  # Calculate the execution time

    print(f"Final Rolling Average Memory Usage: {rolling_average:.3f} MB")
    print(f"Minimum Memory Usage: {min_memory:.3f} MB")
    print(f"Maximum Memory Usage: {max_memory:.3f} MB")
    print(f"Execution Time: {execution_time:.3f} seconds")

def main():
    args = get_args()
    beamsize = str(args.b)
    is_verbose = "1" if args.verbose else "0"
    path = os.path.dirname(os.path.abspath(__file__))
    print(args)

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

    # Monitor memory usage
    monitor_memory_usage(process)

    # Capture and print the output after the process completes
    stdout, stderr = process.communicate()
    if stdout:
        print(stdout.decode())
    if stderr:
        print("There was an error!")
        print(stderr.decode())

if __name__ == "__main__":
    main()
