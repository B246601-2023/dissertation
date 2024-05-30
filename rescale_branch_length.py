import os
from Bio import Phylo
import argparse

def scale_tree_branches_biopython(tree_file, output_file, scale_factor):
    # Load tree
    tree = Phylo.read(tree_file, "newick")

    # Rescale the branch length
    for clade in tree.find_clades():
        if clade.branch_length:
            clade.branch_length *= scale_factor

    # Write the scaled tree to output file
    Phylo.write(tree, output_file, "newick")

def process_trees(input_dir, output_dir, scale_factor):
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Process each file in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".tre"):
            # Generate output filename with scale factor
            base_name = os.path.splitext(filename)[0]  # Remove the extension
            scaled_filename = f"{base_name}_scale{scale_factor}.tre"
            input_file = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, scaled_filename)
            scale_tree_branches_biopython(input_file, output_file, scale_factor)

def main():
    # Setup command line argument parsing
    parser = argparse.ArgumentParser(description="Scale branch lengths of phylogenetic trees in a directory.")
    parser.add_argument('--input_dir', type=str, required=True, help="Directory containing input Newick tree files")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory where output Newick tree files will be saved")
    args = parser.parse_args()

    # Process all tree files in the directory with the given scale factors
    scale_factors = [0.01, 0.1, 0.5, 1, 5, 10, 1000]
    for x in scale_factors:
        process_trees(args.input_dir, args.output_dir, x)

if __name__ == "__main__":
    main()

