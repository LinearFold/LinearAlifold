import argparse
import subprocess
import psutil
import time

def monitor_memory_usage(executable_path, input_file, output_file, args, interval=1):
    start_time = time.time()  # Record the start time
    total_measurements = 0  # Count of measurements
    rolling_average = 0.0  # Rolling average of memory usage
    min_memory = float('inf')  # Initialize min memory to infinity
    max_memory = float('-inf')  # Initialize max memory to negative infinity

    try:
        # Construct the command
        command = [executable_path] + args
        print("Executing command:", ' '.join(command))  # Debugging: print the command

        # Start the executable as a subprocess with stdin, stdout, and stderr redirection
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        # Open the input file and send its contents to the subprocess
        with open(input_file, 'rb') as file:
            process.stdin.write(file.read())
        process.stdin.close()

        # Monitor the memory usage while the process is running
        while True:
            if process.poll() is not None:
                # Process has terminated
                break

            new_memory_usage = psutil.Process(process.pid).memory_info().rss / (1024 * 1024)  # Memory usage in MB
            total_measurements += 1
            rolling_average += (new_memory_usage - rolling_average) / total_measurements
            min_memory = min(min_memory, new_memory_usage)
            max_memory = max(max_memory, new_memory_usage)
            print(f"Current Memory Usage: {new_memory_usage:.3f} MB, Rolling Average: {rolling_average:.3f} MB")
            time.sleep(interval)

        # Capture the output
        with open(output_file, 'wb') as outfile:
            outfile.write(process.stdout.read())

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
    parser = argparse.ArgumentParser(description="Run and monitor an executable.")
    parser.add_argument('-i', '--input', required=True, help='Input file path')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    parser.add_argument('--pt', action='store_true', help='Run Partition mode')

    args = parser.parse_args()

    executable_path = './RNAalifold'

    additional_args = ['-f', 'F', '-r', '--noPS', '--quiet']
    

    if args.pt:
        additional_args += ['-p', '--MEA', '--noDP']

    monitor_memory_usage(executable_path, args.input, args.output, additional_args)

if __name__ == "__main__":
    main()
