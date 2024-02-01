with open('test_data.fasta', 'r') as infile:
    lines = infile.readlines()

output = []

for line in lines:
    line = line.strip()  # remove leading/trailing white space
    if not line:
        continue  # skip blank lines

    if line[0].isdigit():
        # add ">" before sequence identifier
        line = '>' + line
    
    elif line[0] not in ['A', 'C', 'G', 'T', 'U']:
        continue

    output.append(line)
with open('formatted_test_data.txt', 'w') as outfile:
    for line in output[:-1]:
        outfile.write(line + '\n')
