import pyvolve
import os
import argparse

# set a mutation rate
mutation_rates = {
    "AC": 0.039,
    "AG": 0.310,
    "AT": 0.123,
    "CA": 0.140,
    "CG": 0.022,
    "CT": 3.028,
    "GA": 0.747,
    "GC": 0.113,
    "GT": 2.953,
    "TA": 0.056,
    "TC": 0.261,
    "TG": 0.036
}


# function to ajdust branch length
# eg. 0.5x 1x(origin) 10x 100x 
def adjust_branch_lengths(tree_file, factors=[0.01, 0.1, 0.5, 1, 5, 10, 1000]):
    
    scaled_trees = {}
    name, _ = os.path.splitext(os.path.basename(tree_file))
    for factor in factors :
        tree = pyvolve.read_tree(file = tree_file, scale_tree = factor)
        scaled_trees[f"{name}_scale{factor}"] = tree
    
    return scaled_trees

def generate_sequences(tree_dict, out_dir, mutation_rates):
    for name, tree in tree_dict.items():
        model = pyvolve.Model("nucleotide", {"mu": mutation_rates})
        partition = pyvolve.Partition(models=model, size=20000)  # Using a size of 10,000 as an example
        evolver = pyvolve.Evolver(tree=tree, partitions=partition)
        # Ensure output directory exists
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, name)
        evolver(ratefile=f"{out_path}_ratefile.txt", infofile=f"{out_path}_infofile.txt", seqfile=f"{out_path}_seqfile.fasta", write_anc = True)    

def main():
    parser = argparse.ArgumentParser(description="Generate sequences from a phylogenetic tree with specific mutation rates.")
    parser.add_argument('--tree', type=str, required=True, help="Input tree file")
    parser.add_argument('--out', type=str, required=True, help="Output sequence file")
    args = parser.parse_args()

    scaled_trees = adjust_branch_lengths(args.tree)
    generate_sequences(scaled_trees, args.out, mutation_rates)

if __name__ == "__main__":
    main()
