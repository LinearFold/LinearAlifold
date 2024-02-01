import os
import re  # Python's regular expression module
import subprocess
import tempfile
import time

# Make sure the tmp_fasta directory exists
if not os.path.exists('tmp_fasta'):
    os.makedirs('tmp_fasta')

with open("formatted_test_data.txt", "r") as f:
    lines = f.readlines()

total_runtime_lin_alifold = 0
total_runtime_linear_alifold = 0
equal_count = 0
unequal_count = 0

for i in range(0, len(lines), 2):  # Assuming the sequence name and the sequence are always paired
    # Create temp file in the current directory inside tmp_fasta folder
    with tempfile.NamedTemporaryFile(suffix=".fasta", dir='./tmp_fasta', delete=False) as tmp:
        # Add ">" before sequence identifier and include the sequence in the next line
        tmp.write(">".encode() + lines[i].encode())
        seq_name = lines[i].strip()  # save sequence name
        if i+1 < len(lines):
            tmp.write(lines[i+1].encode())
        temp_file_name = tmp.name

    print(f"Processing sequence {seq_name}:")

    # Run the LinAliFold command for each sequence and capture output
    start_time = time.time()
    result_lin_alifold = subprocess.run(['./../tmp/LinAliFold-CentroidLinAliFold/bin/LinAliFold', '-i', temp_file_name], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    runtime_lin_alifold = time.time() - start_time
    total_runtime_lin_alifold += runtime_lin_alifold

    # Extract score from the output using regular expressions
    score_regex = re.compile(r'score:(-?\d+\.\d+)')  # adjust this regex to match the exact format of your scores
    matches_lin_alifold = score_regex.findall(result_lin_alifold.stdout)

    # Run the linear_alifold command for each sequence and capture output
    start_time = time.time()
    result_linear_alifold = subprocess.run(['./../linearalifold_mfe.py', '--em', '2'], check=True, input=open(temp_file_name, 'r').read(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    runtime_linear_alifold = time.time() - start_time
    total_runtime_linear_alifold += runtime_linear_alifold

    # Extract score from the output using regular expressions
    score_regex = re.compile(r'\((-?\d+\.\d+)')
    matches_linear_alifold = score_regex.findall(result_linear_alifold.stdout)

    # Delete temporary file after usage
    os.remove(temp_file_name)

    # Print the score(s)
    for match_lin_alifold, match_linear_alifold in zip(matches_lin_alifold, matches_linear_alifold):
        print("LinAliFold Score: ", match_lin_alifold)
        print("LinearAliFold Score: ", match_linear_alifold)
        if match_lin_alifold == match_linear_alifold:
            equal_count += 1
            print("Scores are equal.")
        else:
            unequal_count += 1
            print("Scores are not equal.")

    # Print runtimes and compute percentage difference
    print(f"LinAliFold runtime: {runtime_lin_alifold} seconds")
    print(f"LinearAliFold runtime: {runtime_linear_alifold} seconds")
    runtime_diff = (runtime_lin_alifold - runtime_linear_alifold) / runtime_linear_alifold * 100
    print(f"Percentage difference in runtime: {runtime_diff}%\n")

print("\nTotal metrics:")
print("Total Sequences Tested:", unequal_count + equal_count)
print("Total runtime of LinAliFold: ", total_runtime_lin_alifold)
print("Total runtime of LinearAliFold: ", total_runtime_linear_alifold)
total_runtime_diff = (total_runtime_lin_alifold - total_runtime_linear_alifold) / total_runtime_linear_alifold * 100
print(f"Total percentage difference in runtime: {total_runtime_diff}%")
print(f"Equal scores count: {equal_count}")
print(f"Unequal scores count: {unequal_count}")
