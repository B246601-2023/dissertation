import pyvolve
import os
import argparse

# set a sample mutation rate
mutation_rates = {
    "AC": 0.5,
    "AG": 2.0,
    "AT": 1.0,
    "CG": 0.7,
    "CT": 1.5,
    "GT": 1.2
}

# function to ajdust mutation rates
# eg. 0.5x 1x(origin) 1.5x 2x 
def adjust_mutation_rates(mutation_rates, factors=[0.5, 1, 1.5, 2]):
    mutation_rate_sets = {}
    for factor in factors:
        adjusted_rates = {k: v * factor for k, v in mutation_rates.items()}
        mutation_rate_sets[f"{factor}x"] = adjusted_rates
    return mutation_rate_sets


def generate_sequences(tree_file, out, mutation_rate_sets):
    
    #read file
    mytree = pyvolve.read_tree(file = tree_file)
    
    freqs = [0.25, 0.25, 0.25, 0.25]
    
    for name, rates in mutation_rate_sets.items():
        # set a model
        mymodel = pyvolve.Model(
            "nucleotide",
            {"mu": rates, "state_freqs": freqs},
            alpha=0.5, num_categories=3, pinv=0.25
        )

        my_partition = pyvolve.Partition(models=mymodel, size=100) # save time, should be 20000
        evolver = pyvolve.Evolver(tree=mytree, partitions=my_partition)
        

        os.makedirs(os.path.dirname(out), exist_ok=True)

        # generate file and save
        evolver(ratefile=f"{out}_{name}_ratefile.txt", infofile=f"{out}_{name}_infofile.txt", seqfile=f"{out}_{name}_seqfile.fasta")

    

def main():
    parser = argparse.ArgumentParser(description="Generate sequences from a phylogenetic tree with specific mutation rates.")
    parser.add_argument('--tree', type=str, required=True, help="Input tree file")
    parser.add_argument('--out', type=str, required=True, help="Output sequence file")
    args = parser.parse_args()

    mutation_rate_sets = adjust_mutation_rates(mutation_rates)

    generate_sequences(args.tree, args.out, mutation_rate_sets)

if __name__ == "__main__":
    main()
