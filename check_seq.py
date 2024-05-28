from Bio import SeqIO
import os
import dendropy
import pandas as pd
import argparse
import matplotlib.pyplot as plt

def calculate_phylogenetic_distance(input_path):
    data = []
    tree = dendropy.Tree.get(
        path=input_path,
        schema="newick",
        tree_offset=None,
        )
    pdc = tree.phylogenetic_distance_matrix()
    for i, t1 in enumerate(tree.taxon_namespace[:-1]):
        for t2 in tree.taxon_namespace[i+1:]:
            data.append({
                'Taxon1': t1.label,
                'Taxon2': t2.label,
                'Phylogenetic Distance': pdc(t1, t2)
            })
    
    df = pd.DataFrame(data)

    # 根据'Phylogenetic Distance'排序数据
    df_sorted = df.sort_values(by='Phylogenetic Distance')

    # 等距选择20个样本
    if len(df_sorted) > 20:
        step = len(df_sorted) // 20
        df_sampled = df_sorted.iloc[::step].head(20)
    else:
        df_sampled = df_sorted
    return df_sampled

def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of the same length to compute Hamming distance")
    
    distance = sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))
    return distance


def load_sequences(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def calculate_distances(fasta_file, df):

    sequences = load_sequences(fasta_file)

    # calculate hamming distance
    results = []
    for _, row in df.iterrows():
        taxon1 = row['Taxon1']
        taxon2 = row['Taxon2']
        if taxon1 in sequences and taxon2 in sequences:
            distance = hamming_distance(sequences[taxon1], sequences[taxon2])
            results.append({
                'Taxon1': taxon1,
                'Taxon2': taxon2,
                'Hamming Distance': distance,
                'Phylogenetic Distance': row['Phylogenetic Distance']
            })
        else:
            results.append({
                'Taxon1': taxon1,
                'Taxon2': taxon2,
                'Hamming Distance': None,
                'Phylogenetic Distance': row['Phylogenetic Distance']
            })

    return pd.DataFrame(results)

def plot_distances(input,df):
    name, _ = os.path.splitext(os.path.basename(input))
    df = df.dropna()  # Remove any rows with None values
    plt.scatter(df['Phylogenetic Distance'], df['Hamming Distance'])
    plt.xlabel('Phylogenetic Distance')
    plt.ylabel('Hamming Distance')
    plt.title('Phylogenetic vs Hamming Distance')
    plt.grid(True)
    plt.savefig(f"{name}_phplot.pdf")

def main():
    parser = argparse.ArgumentParser(description="Calculate phylogenetic distances between all pairs of taxa in a tree and save to CSV.")
    parser.add_argument('--tree', type=str, help="Input path of the Newick tree file")
    parser.add_argument('--fasta', type=str, help="Input path of the fasta file")

    #parser.add_argument('output_path', type=str, help="Output path for the CSV file")
    args = parser.parse_args()
    df = calculate_phylogenetic_distance(args.tree)
    res = calculate_distances(args.fasta,df)
    plot_distances(args.tree,res)

if __name__ == "__main__":
    main()