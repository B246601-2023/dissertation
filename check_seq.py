from Bio import SeqIO
import os
import dendropy
import pandas as pd
import argparse
import matplotlib.pyplot as plt

def getTree(input_path):
    tree = dendropy.Tree.get(
        path=input_path,
        schema="newick",
        tree_offset=None,
        )
    return tree

def calculate_phylogenetic_distance(tree):
    data = []
    pdc = tree.phylogenetic_distance_matrix()
    for i, t1 in enumerate(tree.taxon_namespace[:-1]):
        for t2 in tree.taxon_namespace[i+1:]:
            data.append({
                'Taxon1': t1.label,
                'Taxon2': t2.label,
                'Phylogenetic Distance': pdc(t1, t2)
            })
    
    df = pd.DataFrame(data)
    df_sorted = df.sort_values(by='Phylogenetic Distance')
    
    #choose 20 sample for Phylogenetic Distance-Hamming plot
    if len(df_sorted) > 20:
        step = len(df_sorted) // 20
        df_sampled = df_sorted.iloc[::step].head(20)
    else:
        df_sampled = df_sorted
    
    # 20 samples
    df_box = df.sort_values(by='Phylogenetic Distance', ascending=False).head(20)
    return df_sampled,df_box

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

def plot_distances(name, df_sample, df_box, scales,fasta_dir):
    fig, ax = plt.subplots(1, 2, figsize=(14, 7))
    
    hamming = calculate_distances(f"{fasta_dir}/{name}_scale1.0_seqfile.fasta", df_sample)
    # Plot phylogenetic vs Hamming distance
    ax[0].scatter(hamming['Phylogenetic Distance'], hamming['Hamming Distance'])
    ax[0].set_xscale("log")
    ax[0].set_xlabel('Phylogenetic Distance')
    ax[0].set_ylabel('Hamming Distance')
    ax[0].set_title('Phylogenetic vs Hamming Distance')
    ax[0].grid(True)

    # Boxplot of different scales
    distance = []
    for x in scales :
        df = calculate_distances(f"{fasta_dir}/{name}_scale{x}_seqfile.fasta", df_box)
        distance.append([d for d in df['Hamming Distance'] if d is not None])
    
    scale_labels = [f"{x}x" for x in scales]
    #distances = [calculate_distances(f"{fasta_dir}/{name}_scale{x}_seqfile.fasta", df_box) for x in scales]
    ax[1].boxplot(distance, patch_artist=True, labels=scale_labels)
    ax[1].set_xlabel('Scale Factor')
    ax[1].set_ylabel('Hamming Distance')
    ax[1].set_title('Hamming Distance at Different Scales')
    ax[1].grid(True)

    plt.savefig(f"{name}_phplot.pdf")

def main():
    parser = argparse.ArgumentParser(description="Calculate phylogenetic distances between all pairs of taxa in a tree and save to CSV.")
    parser.add_argument('--tree', type=str, help="Input path of the Newick tree file")
    parser.add_argument('--fasta_dir', type=str, help="Input path of the fasta file")

    #parser.add_argument('output_path', type=str, help="Output path for the CSV file")
    args = parser.parse_args()
    tree = getTree(args.tree)
    df_sample,df_box = calculate_phylogenetic_distance(tree)
    scales = [0.0001, 0.001, 0.01, 0.1, 1.0, 10.0]
    name, _ = os.path.splitext(os.path.basename(args.tree))
    plot_distances(name,df_sample,df_box,scales,args.fasta_dir)

if __name__ == "__main__":
    main()