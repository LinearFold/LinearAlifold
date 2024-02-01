import argparse
import subprocess
import psutil
import time

def monitor_cpp_memory_usage(executable_path, args, output_file, interval=1):
    start_time = time.time()  # Record the start time
    total_measurements = 0  # Count of measurements
    rolling_average = 0.0  # Rolling average of memory usage
    min_memory = float('inf')  # Initialize min memory to infinity
    max_memory = float('-inf')  # Initialize max memory to negative infinity

    try:
        # Construct the command
        command = [executable_path] + args
        print("Executing command:", ' '.join(command))  # Debugging: print the command

        # Open the output file
        with open(output_file, 'w') as outfile:
            # Start the C++ executable as a subprocess with stdout and stderr redirected
            process = subprocess.Popen(command, stdout=outfile, stderr=subprocess.STDOUT)

            # Wait for the process to start
            # time.sleep(1)

            # Check if the process has terminated early
            if process.poll() is not None:
                print("Process terminated early. Return code:", process.returncode)
                return

            # Get the process using psutil
            cpp_process = psutil.Process(process.pid)

            # Monitor the memory usage
            while True:
                # Check if the process has terminated
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
    except KeyboardInterrupt:
        print("Monitoring interrupted")
    except Exception as e:
        print(f"An error occurred: {e}")

    end_time = time.time()  # Record the end time
    execution_time = end_time - start_time  # Calculate the execution time

    print(f"Final Rolling Average Memory Usage: {rolling_average:.3f} MB")
    print(f"Minimum Memory Usage: {min_memory:.3f} MB")
    print(f"Maximum Memory Usage: {max_memory:.3f} MB")
    print(f"Execution Time: {execution_time:.3f} seconds")

def main():
    parser = argparse.ArgumentParser(description="Run and monitor C++ executable.")
    parser.add_argument('-i', '--input', required=True, help='Input file path')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    parser.add_argument('--cen', action='store_true', help='Run CentroidLinAliFold')

    args = parser.parse_args()

    if args.cen:
        executable_path = './bin/CentroidLinAliFold'
        cpp_args = ['-i', args.input, '-p', '1', '-o', '1', '-r', '1']
    else:
        executable_path = './bin/LinAliFold'
        cpp_args = ['-i', args.input, '-r', '1']

    monitor_cpp_memory_usage(executable_path, cpp_args, args.output)

if __name__ == "__main__":
    main()



