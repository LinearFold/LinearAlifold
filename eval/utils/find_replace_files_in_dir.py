import os
import argparse

def replace_string_in_files(directory, old_string, new_string):
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path):
            with open(file_path, 'r') as file:
                file_contents = file.read()

            file_contents = file_contents.replace(old_string, new_string)

            with open(file_path, 'w') as file:
                file.write(file_contents)
            print(f"Updated file: {file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Replace a specified string in all files of a directory.")
    parser.add_argument("directory", help="Path to the directory containing the files")
    parser.add_argument("old_string", help="The string to be replaced")
    parser.add_argument("new_string", help="The new string to replace the old one", default="")

    args = parser.parse_args()

    replace_string_in_files(args.directory, args.old_string, args.new_string)
