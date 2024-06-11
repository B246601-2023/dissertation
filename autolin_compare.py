import pandas as pd
import argparse
import os
from Bio import Phylo

def handle_data(input_file, output_file):
    """
    Reads a CSV file, processes it by removing 'auto.' prefixes, and saves it to another file.
    
    Parameters:
        input_file (str): Path to the input CSV file.
        output_file (str): Path to the output CSV file where the processed data will be saved.
    """
    columns = [0, 1]
    df = pd.read_csv(input_file, sep='\t', usecols=columns)
    df.columns = ['clade', 'root_id']  # Assume column 0 is old names, column 1 is new names
    df['clade'] = df['clade'].apply(lambda x: x.replace("auto.", ""))
    df.to_csv(output_file, index=False)
    return df

def rename_node(tree, name_pairs):
    """
    Modifies the names of nodes in a phylogenetic tree based on a DataFrame of name pairs.
    
    Parameters:
        tree: Phylo tree object
        name_pairs: DataFrame containing old and new names
    """
    for _, row in name_pairs.iterrows():
        for clade in tree.find_clades():
            if clade.name == row['root_id']:
                clade.name = row['clade']
    return tree

def main():
    parser = argparse.ArgumentParser(description='Process and save data.')
    parser.add_argument('--input', type=str, help='Input file path for auto annotations')
    parser.add_argument('--output', type=str, help='Output file path for modified auto annotations')
    parser.add_argument('--tree_dir', type=str, help='Input tree directory path for annotations')
    args = parser.parse_args()

    # Process the annotation file to get name pairs
    name_df = handle_data(args.input, args.output)
    
    # Assuming the tree name is derived from the input file name
    name, _ = os.path.splitext(os.path.basename(args.input))
    tree_name = name + ".nh"
    tree_path = os.path.join(args.tree_dir, tree_name)
    
    # Read the tree, apply renaming, and save the modified tree
    tree = Phylo.read(tree_path, 'newick')
    new_tree = rename_node(tree, name_df)
    
    # Save the modified tree
    os.makedirs("modified_tree", exist_ok=True)
    modified_tree_path = os.path.join("modified_tree", tree_name)
    Phylo.write(new_tree, modified_tree_path, 'newick')

if __name__ == "__main__":
    main()
