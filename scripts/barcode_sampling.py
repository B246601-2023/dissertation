import sys
import os
import random
import argparse
import pandas as pd
import numpy as np


def main(args):

    # make output dir if needed
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)
        print(f"Random seed set to: {args.seed}")
    else:
        np.random.seed(None)

    # Read barcodes file
    try:
        bcs = pd.read_csv(args.barcodes)
        print(f"Barcode dataset loaded from: {args.barcodes}")
    except FileNotFoundError:
        print(f"Error: The file {args.barcodes} does not exist.")
        sys.exit(1)

    for i in range(1, args.replicates+1):
        # If not fixed num_lineages, sample num_lineages from uniform dist
        if not args.num_lineages:
            n = np.random.randint(args.min_lineages, args.max_lineages + 1)
        else:
            n = args.num_lineages

        print(f"Replicate {i} sampled {n} lineages...")

        # Sample lineages and assign frequencies
        lin_freqs = sample_lineages(bcs, n, args.dist,
                                    args.alpha)

        # Compute SNV frequencies
        snv_freqs = compute_snv_frequencies(bcs, lin_freqs)

        # Convert SNV frequencies to pandas DataFrame
        snv_df = to_variants_table(snv_freqs, reference=args.ref,
                                   fixed_depth=args.depth)

        # Create output file names
        variant_file = f"{args.out}_rep{i+1}.variants.tsv"
        lineage_file = f"{args.out}_rep{i+1}.known_lineages.tsv"

        # Write the variants table to file
        snv_df.to_csv(variant_file, sep='\t', index=False)

        # Write the known lineage frequencies to file
        lin_freqs_df = pd.DataFrame.from_dict(lin_freqs, orient='index',
                                              columns=['Frequency'])
        lin_freqs_df.index.name = 'Lineage'
        lin_freqs_df.to_csv(lineage_file, sep='\t')

    print(f"Done! Wrote outputs using prefix {args.out}")


def to_variants_table(snv_freqs, reference="SPOOF", fixed_depth=1000):
    """
    Converts the SNV frequencies dictionary to a variants table input for 
    Freyja.

    Parameters:
    snv_freqs (dict): A dictionary containing SNV positions and their
                      frequencies.
    reference (str): The reference genome name. Default is "SPOOF".
    fixed_depth (int): The total depth used to calculate DP values.
                       Default is 1000.

    Returns:
    pd.DataFrame: A DataFrame containing the SNV frequencies formatted for
                  Freyja.
    """
    data = []
    for pos, freqs in snv_freqs.items():
        ref = freqs['REF']
        total_alt_freq = sum(freqs['ALT'].values())
        ref_dp = fixed_depth * (1 - total_alt_freq)
        for alt, alt_freq in freqs['ALT'].items():
            alt_dp = fixed_depth * alt_freq
            data.append({
                'REGION': reference,
                'POS': pos,
                'REF': ref,
                'ALT': alt,
                'REF_DP': ref_dp,
                'REF_RV': 0,  # Spoof
                'REF_QUAL': 33,  # Spoof
                'ALT_DP': alt_dp,
                'ALT_RV': 0,  # Spoof
                'ALT_QUAL': 33,  # Spoof
                'ALT_FREQ': alt_freq,
                'TOTAL_DP': fixed_depth,
                'PVAL': 0,  # Spoof
                'PASS': 'TRUE',  # Spoof
                'GFF_FEATURE': 'NA',  # Spoof
                'REF_CODON': 'NA',  # Spoof
                'REF_AA': 'NA',  # Spoof
                'ALT_CODON': 'NA',  # Spoof
                'ALT_AA': 'NA',  # Spoof
            })

    snv_df = pd.DataFrame(data)
    snv_df = snv_df.sort_values(by='POS')

    return snv_df


def compute_snv_frequencies(bcs, lineage_freqs, filter=True):
    """
    Computes the SNV frequencies based on the sampled lineage frequencies.

    Parameters:
    bcs (pd.DataFrame): The barcodes dataframe.
    lineage_freqs (dict): A dictionary containing sampled lineage names and
                          their frequencies.
    filter (bool): Whether to filter out positions where the sum of ALT
                   frequencies is 0. Default is True.

    Returns:
    dict: A dictionary containing SNV positions and their frequencies.
    """
    snv_freqs = {}

    for lineage, freq in lineage_freqs.items():
        lineage_snv = bcs[bcs.iloc[:, 0] == lineage].iloc[:, 1:].to_dict(
            orient='records')[0]
        for snv, snv_freq in lineage_snv.items():
            pos = int(snv[1:-1])
            ref = snv[0]
            alt = snv[-1]
            if pos not in snv_freqs:
                snv_freqs[pos] = {'REF': ref, 'ALT': {}}
            else:
                if snv_freqs[pos]['REF'] != ref:
                    print(f"Warning: Inconsistent REF state at position {pos}"
                          f"Expected {snv_freqs[pos]['REF']}, found {ref}")

            if alt not in snv_freqs[pos]['ALT']:
                snv_freqs[pos]['ALT'][alt] = 0
            snv_freqs[pos]['ALT'][alt] += snv_freq * freq

    # if filter is True, remove sites where only state is ref
    if filter:
        snv_freqs = {
            pos: data
            for pos, data in snv_freqs.items()
            if sum(data['ALT'].values()) > 0
        }
    # Sort the dictionary by position
    snv_freqs = dict(sorted(snv_freqs.items()))

    return snv_freqs


def sample_lineages(bcs, num_lineages, frequency_distribution, alpha):
    """
    Samples lineages from the given barcodes dataframe and assigns frequencies
    based on the specified distribution type.

    Parameters:
    bcs (pd.DataFrame): The barcodes dataframe.
    num_lineages (int): The number of lineages to sample.
    frequency_distribution (str): The type of frequency distribution to use
                                  ('uniform' or 'dirichlet').
    alpha (float): The alpha parameter for the Dirichlet distribution
                   when 'dirichlet' is chosen.

    Returns:
    dict: A dictionary containing the sampled lineage names and the assigned
          frequencies.

    Raises:
    ValueError: If an invalid frequency distribution type is provided.
    """
    sampled_bcs = bcs.sample(num_lineages)

    if frequency_distribution == 'uniform':
        frequencies = np.ones(num_lineages) / num_lineages  # Balanced
    elif frequency_distribution == 'dirichlet':
        alpha = np.linspace(alpha, 1, num_lineages)  # Skewed
        frequencies = np.random.dirichlet(alpha)
    else:
        raise ValueError(
            f"Invalid frequency distribution type: {frequency_distribution}")

    freqs = dict(zip(sampled_bcs.iloc[:, 0], frequencies))

    return freqs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sample from Freyja barcode database file")

    parser.add_argument(
        '-s', '--seed', type=int,
        help='Random seed for reproducibility'
    )
    parser.add_argument(
        '-b', '--barcodes', type=str, required=True,
        help='Freyja barcodes database file'
    )
    parser.add_argument(
        '-n', '--num_lineages', type=int, required=False,
        help='Fixed number of lineages to sample'
    )
    parser.add_argument(
        '-m', '--min_lineages', type=int, default=2,
        help='Minimum number of lineages to sample (for range sampling)'
    )
    parser.add_argument(
        '-M', '--max_lineages', type=int, default=10,
        help='Maximum number of lineages to sample (for range sampling)'
    )
    parser.add_argument(
        '-f', '--dist', choices=['uniform', 'skewed'],
        default='uniform',
        help='Distribution of frequencies assigned to sampled lineages'
    )
    parser.add_argument(
        '-a', '--alpha', type=float, default=2.0,
        help='Alpha parameter for Dirichlet distribution'
    )
    parser.add_argument(
        '-r', '--replicates', type=int, default=1,
        help='Number of replicates for sampling'
    )
    parser.add_argument(
        '-o', '--out', type=str, default="output/out",
        help='Output file prefix'
    )
    parser.add_argument(
        '-R', '--ref', type=str, default="NC_045512.2",
        help="Reference genome name for output files")
    parser.add_argument(
        '-d', '--depth', type=int, default=1000,
        help="Sequencing depth to use for output files")

    args = parser.parse_args()
    main(args)