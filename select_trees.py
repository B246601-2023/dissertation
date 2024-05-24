import os
import re
import random
import argparse

# Currently, trees are chosen randomly

def select_and_clean_trees(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # List all tree files
    tree_files = [f for f in os.listdir(input_dir) if f.endswith('.tre')]

    if not tree_files:
        print("No tree files found in the input directory.")
        return

    # Extract unique SD values
    sd_values = set(re.findall(r'_s\d+\.\d+', ' '.join(tree_files)))

    for sd in sd_values:
        # Select a random tree file for each SD value
        sd_tree_files = [f for f in tree_files if sd in f]
        selected_tree = random.choice(sd_tree_files)

        # Clean the tree file
        with open(os.path.join(input_dir, selected_tree), 'r') as infile:
            content = infile.read()
        cleaned_content = content.replace('[&R] ', '')

        # Save the cleaned tree file to the output directory
        with open(os.path.join(output_dir,selected_tree), 'w') as outfile:
            outfile.write(cleaned_content)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Select and clean tree files.")
    parser.add_argument('--input', type=str, required=True, help="Input directory containing tree files")
    parser.add_argument('--out', type=str, required=True, help="Output directory for cleaned tree files")
    args = parser.parse_args()

    select_and_clean_trees(args.input, args.out)
