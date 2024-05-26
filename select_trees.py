import os
import re
import random
import argparse

# Currently, trees are chosen randomly

def select_and_clean_trees(input_dir, output_dir, selection_file):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    try:
        with open(os.path.join(input_dir, selection_file), 'r') as file:
            selected_tree_files = [line.strip() for line in file if line.strip()]
    except FileNotFoundError:
        print(f"The file {selection_file} was not found in the input directory.")
        return

    for filename_long in selected_tree_files:
        # read each selected tree from input dir
        filename = os.path.basename(filename_long)
        file_path = os.path.join(input_dir, filename)
        if not os.path.exists(file_path):
            print(file_path)
            print(f"File {filename} listed in {selection_file} does not exist in the input directory.")
            continue

        # read and clean trees
        with open(file_path, 'r') as infile:
            content = infile.read()
        cleaned_content = content.replace('[&R] ', '')

        # Save the cleaned tree file to the output directory
        with open(os.path.join(output_dir,filename), 'w') as outfile:
            outfile.write(cleaned_content)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Select and clean tree files.")
    parser.add_argument('--input', type=str, required=True, help="Input directory containing tree files")
    parser.add_argument('--out', type=str, required=True, help="Output directory for cleaned tree files")
    parser.add_argument('--select', type=str, default='trees_selected_trees', help="File listing the trees to be processed")
    args = parser.parse_args()

    select_and_clean_trees(args.input, args.out, args.select)
