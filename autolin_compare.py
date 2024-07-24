import pandas as pd
import argparse
import os
from Bio import Phylo

def handle_data(sample, clade, output_file):
    """
    Reads a CSV file, processes it by removing 'auto.' prefixes, and saves it to another file.
    
    Parameters:
        sample (str): Path to the input sample path file.
        clade (str): Path to the input clade path file.
        output_file (str): Path to the output CSV file where the processed data will be saved.
    """
    columns = [0, 1]
    df_sample = pd.read_csv(sample, sep='\t', header=None, names=["root_id","path"])
    df_clade = pd.read_csv(clade, sep="\t")

    df_clade_clean = df_clade[df_clade['root_id'].isin(df_sample['root_id'])]
    
    df_clade['clade'] = df_clade['clade'].apply(lambda x: x.replace("auto.", "") if "auto." in x else x)
    
    df_clade_clean.to_csv(output_file, sep="\t", index=False)

    return df_clade

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
    parser.add_argument('--sample', type=str, help='Input file path for auto sample annotations')
    parser.add_argument('--clade', type=str, help='Input file path for auto clade annotations')
    parser.add_argument('--output', type=str, help='Output file path for modified auto annotations')
    parser.add_argument('--tree', type=str, help='Input tree directory path for annotations')
    args = parser.parse_args()

    # Process the annotation file to get name pairs
    name_df = handle_data(args.sample,args.clade, args.output)
    
    # Assuming the tree name is derived from the input file name
    # name, _ = os.path.splitext(os.path.basename(args.input))
    # tree_name = name + ".nh"
    # tree_path = os.path.join(args.tree_dir, tree_name)
    
    # Read the tree, apply renaming, and save the modified tree
    tree = Phylo.read(args.tree, 'newick')
    new_tree = rename_node(tree, name_df)
    
    # Save the modified tree
    os.makedirs("modified_tree", exist_ok=True)
    modified_tree_path = os.path.join("modified_tree", args.tree)
    Phylo.write(new_tree, modified_tree_path, 'newick')

if __name__ == "__main__":
    main()
