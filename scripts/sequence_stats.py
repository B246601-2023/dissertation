import os
import pandas as pd
from Bio import AlignIO
import argparse
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import subprocess
import tempfile
import struct


def run_dashing2_cmp(sketch_file, k=31):
    output_dist = sketch_file + ".dist"
    
    try:
        result = subprocess.run(
            ['dashing2', 'cmp',  
             '-k', str(k), 
             '--asymmetric-all-pairs', 
             '--compute-edit-distance',
             sketch_file, 
             '-o', output_dist],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running dashing2 cmp: {e}")
        print("Command output (cmp):")
        print("STDOUT:", e.stdout.decode())
        print("STDERR:", e.stderr.decode())
        raise
    
    return output_dist

def parse_dashing_output(output_dist):
    dist_matrix = []
    try:
        with open(output_dist, 'rb') as f:
            while True:
                chunk = f.read(4)
                if not chunk:
                    break
                num = struct.unpack('f', chunk)[0]
                dist_matrix.append(num)
    except Exception as e:
        print(f"Error reading the distance matrix: {e}")
        raise
    num_elements = len(dist_matrix)
    size = int((1 + (1 + 8 * num_elements)**0.5) / 2)
    square_matrix = np.zeros((size, size), dtype=float)
    idx = 0
    for i in range(size):
        for j in range(i + 1, size):
            square_matrix[i, j] = dist_matrix[idx]
            square_matrix[j, i] = dist_matrix[idx]
            idx += 1
    return square_matrix

def get_pairwise_stats(distance_matrix):
    distances = distance_matrix[np.tril_indices_from(distance_matrix, k=-1)]
    mean_distance = np.nanmean(distances)
    median_distance = np.nanmedian(distances)
    non_zero_distances = distances[distances > 0]
    min_distance = np.min(non_zero_distances) if len(non_zero_distances) > 0 else np.nan

    if len(non_zero_distances) > 0:
        min_idx = np.unravel_index(np.nanargmin(np.where(distance_matrix > 0, distance_matrix, np.nan)), distance_matrix.shape)
        min_pair = (min_idx[0], min_idx[1])
    else:
        min_pair = (0, 0)

    return mean_distance, median_distance, min_distance, (int(min_pair[0]), int(min_pair[1]))

def calculate_p_distance(seq1, seq2):
    hamming_distance = calculate_hamming_distance(seq1, seq2)
    seq_length = max(len(seq1), len(seq2))
    return hamming_distance / seq_length

def calculate_hamming_distance(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def calc_haplotype_diversity(alignment):
    sequences = [str(record.seq) for record in alignment]
    unique_sequences = set(sequences)
    n = len(sequences)
    k = len(unique_sequences)
    haplotype_frequencies = [sequences.count(haplotype) / n for haplotype in unique_sequences]
    haplotype_diversity = 1 - sum(frequency ** 2 for frequency in haplotype_frequencies)
    return haplotype_diversity

def calculate_snp_stats(directory, k=31):
    data = []

    for filename in os.listdir(directory):
        if filename.endswith('.fasta'):
            parts = filename.split('_')
            SD = float(parts[1][1:])
            REP = int(parts[2][1:])
            SCALE = float(parts[3][5:])
            filepath = os.path.join(directory, filename)
            alignment = AlignIO.read(filepath, 'fasta')
            snp_count = 0
            alignment_length = alignment.get_alignment_length()
            for i in range(alignment_length):
                column = alignment[:, i]
                if len(set(column)) > 1:
                    snp_count += 1
            proportion_variable = snp_count / alignment_length if alignment_length > 0 else 0
            haplotype_diversity = calc_haplotype_diversity(alignment)
            sequences = [str(record.seq) for record in alignment]
            sequence_counts = Counter(sequences)
            num_unique_haplotypes = len(sequence_counts)
            num_duplicates = sum(count - 1 for count in sequence_counts.values() if count > 1)
            
            # Run Dashing2 for sketching and distance calculation
            dist_file = run_dashing2_cmp(filepath, k)
            distance_matrix = parse_dashing_output(dist_file)
            
            mean_kmer_distance, median_kmer_distance, min_kmer_distance, min_pair = get_pairwise_stats(distance_matrix)
            min_p_distance = calculate_p_distance(alignment[min_pair[0]].seq, alignment[min_pair[1]].seq)
            min_hamming = calculate_hamming_distance(alignment[min_pair[0]].seq, alignment[min_pair[1]].seq)

            data.append({
                'SD': SD, 
                'REP': REP, 
                'SCALE': SCALE, 
                'SNPs': snp_count, 
                'Num_Sequences': len(alignment),
                'Length': alignment_length,
                'Prop_Variable': proportion_variable,
                'Haplotype_Diversity': haplotype_diversity,
                'Num_Duplicate_Haplotypes': num_duplicates,
                'Num_Unique_Haplotypes': num_unique_haplotypes,
                'Median_Edit_Dist': median_kmer_distance,
                'Min_Nonzero_P_Distance': min_p_distance,
                'Min_Nonzero_Hamming_Distance': min_hamming,
            })

    df = pd.DataFrame(data)
    df.to_csv('snp_stats.tsv', index=False, sep="\t")
    return df

def plot_stats(df):
    variables = ['SNPs', 'Prop_Variable', 'Haplotype_Diversity', 
                 'Num_Duplicate_Haplotypes', 'Num_Unique_Haplotypes', 
                 'Median_Edit_Dist', 'Min_Nonzero_P_Distance',
                 'Min_Nonzero_Hamming_Distance']

    scales = df['SCALE'].unique()
    num_vars = len(variables)

    fig, axes = plt.subplots(nrows=len(scales), ncols=num_vars, figsize=(num_vars * 5, len(scales) * 5))
    for i, scale in enumerate(scales):
        for j, variable in enumerate(variables):
            ax = axes[i, j]
            sns.boxplot(data=df[df['SCALE'] == scale], x='SD', y=variable, hue='REP', ax=ax)
            ax.set_title(f'SCALE = {scale}, Variable = {variable}')
            ax.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig('snp_stats.pdf')

def main():
    parser = argparse.ArgumentParser(description="Calculate SNP statistics from fasta files")
    parser.add_argument('--directory', type=str, default='results/sequences', help='Directory containing fasta files')
    parser.add_argument('--kmer_len', type=int, default=31, help='Kmer len for sketching')
    args = parser.parse_args()

    df = calculate_snp_stats(args.directory, args.kmer_len)

    plot_stats(df)

if __name__ == '__main__':
    main()