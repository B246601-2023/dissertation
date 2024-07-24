import sys
import os
import random
import argparse
import pandas as pd
import numpy as np
import math


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
    # try:
    #     bcs = pd.read_csv(args.barcodes)
    #     if bcs.columns[0] == "Unnamed: 0.1" :
    #         bcs = bcs.drop(bcs.columns[0], axis=1)
    #         bcs.rename(columns={'Unnamed: 0': 'index'}, inplace=True)
    #     print(f"Barcode dataset loaded from: {args.barcodes}")
    # except FileNotFoundError:
    #     print(f"Error: The file {args.barcodes} does not exist.")
    #     sys.exit(1)
    bcs, bcs_recombine = read_and_process_barcodes(args.barcodes)

    for i in range(1, args.replicates+1):
        # If not fixed num_lineages, sample num_lineages from uniform dist
        if not args.num_lineages:
            n = np.random.randint(args.min_lineages, args.max_lineages + 1)
        else:
            n = args.num_lineages

        print(f"Replicate {i} sampled {n} lineages...")

        # Sample lineages and assign frequencies
        if args.high :
            bcs = pd.concat([bcs, bcs_recombine])
            lin_freqs = sample_lineages(bcs, n, args.dist,
                                    args.alpha, False)
        elif args.low :
            lineage_num = bcs['index'].nunique()
            recombinant_num_low = math.ceil((0.12 * lineage_num)/(1-0.12))
            bcs_low = bcs_recombine.sample(n=recombinant_num_low, random_state=3)
            bcs = pd.concat([bcs, bcs_low])
            lin_freqs = sample_lineages(bcs, n, args.dist, args.alpha, False)
        else :
            lin_freqs = sample_lineages(bcs, n, args.dist,
                                    args.alpha, args.cluster)

        # Compute SNV frequencies
        try:
            snv_freqs = compute_snv_frequencies(bcs, lin_freqs,False)
        except IndexError as e:
            variant_file = f"empty.variants.tsv"
            lineage_file = f"empty.known_lineages.tsv"
            with open(variant_file,'w') as file :
                file.write("empty")
            lin_freqs_df = pd.DataFrame.from_dict(lin_freqs, orient='index',
                                              columns=['Frequency'])
            lin_freqs_df.index.name = 'Lineage'
            lin_freqs_df.to_csv(lineage_file, sep='\t')
            continue

        # Convert SNV frequencies to pandas DataFrame
        snv_df = to_variants_table(snv_freqs, reference=args.ref,
                                   fixed_depth=args.depth)

        # Create output file names
        variant_file = f"{args.out}_rep{i}.variants.tsv"
        lineage_file = f"{args.out}_rep{i}.known_lineages.tsv"

        # Write the variants table to file
        snv_df.to_csv(variant_file, sep='\t', index=False)

        # Write the known lineage frequencies to file
        lin_freqs_df = pd.DataFrame.from_dict(lin_freqs, orient='index',
                                              columns=['Frequency'])
        lin_freqs_df.index.name = 'Lineage'
        lin_freqs_df.to_csv(lineage_file, sep='\t')

    print(f"Done! Wrote outputs using prefix {args.out}")


def read_and_process_barcodes(filepath):
    """
    Reads a barcode file, handles unexpected index columns, and separates rows
    with 'X' in the 'index' column.

    Parameters:
    filepath (str): Path to the barcode CSV file.

    Returns:
    tuple: A tuple containing the cleaned main DataFrame and a DataFrame of rows
           with 'X' in the 'index' column.
    """
    try:
        # read files
        bcs = pd.read_csv(filepath)
        print(f"Barcode dataset loaded from: {filepath}")

        # check if 'unnamed' exsists and deal it
        if bcs.columns[0] == "Unnamed: 0.1":
            bcs.drop(bcs.columns[0], axis=1, inplace=True)
        if 'Unnamed: 0' in bcs.columns:
            bcs.rename(columns={'Unnamed: 0': 'index'}, inplace=True)

        # delete row with 'X'
        rows_with_x = bcs[bcs['index'].str.contains('X', na=False)]
        bcs = bcs[~bcs['index'].str.contains('X', na=False)]

        return bcs, rows_with_x

    except FileNotFoundError:
        print(f"Error: The file {filepath} does not exist.")
        sys.exit(1)

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
        # Check if 'POS' column exists before sorting
    if 'POS' in snv_df.columns:
        snv_df = snv_df.sort_values(by='POS')
    else:
        print("Error: 'POS' column not found in the DataFrame.")
    #snv_df = snv_df.sort_values(by='POS')

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
        # # exit when barcode is null
        # filtered_bcs = bcs[bcs.iloc[:, 0] == lineage]
        # if filtered_bcs.empty or filtered_bcs.columns.size == 0:
        #     print(f"Error: No data found for lineage: {lineage}")
        #     sys.exit(0)

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
    #print("snv_freqs:")
    #print(snv_freqs)
    return snv_freqs


def sample_lineages(bcs, num_lineages, frequency_distribution, alpha, cluster):
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
    if cluster :
        bcs_clean = bcs[~bcs['index'].str.contains('X', na=False)]
        bcs_clean['group'] = bcs_clean['index'].apply(lambda x: x.split('_')[0])
        grouped_bcs = {k: v.drop('group', axis=1) for k, v in bcs_clean.groupby('group')}
        bcs_sets = [df for df in grouped_bcs.values() if len(df) > num_lineages]
        if bcs_sets:
            sampled_bcs = random.choice(bcs_sets)
        else:
            bcs_sets_2 = [df for df in grouped_bcs.values() if len(df) > 1]
            sampled_bcs = random.choice(bcs_sets_2)
            num_lineages = len(sampled_bcs)

    else :
        sampled_bcs = bcs.sample(num_lineages)

    if frequency_distribution == 'uniform':
        frequencies = np.ones(num_lineages) / num_lineages  # Balanced
    elif frequency_distribution == 'skewed':
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
        '-H', '--high',
        action='store_true',  
        help='Enable high rate recombinants in samples'
    )
    parser.add_argument(
        '-l', '--low',
        action='store_true',  
    help='Enable low rate recombinants in samples'
    )
    parser.add_argument(
        '-c', '--cluster',
        action='store_true',  
        help='Enable phylogenetically clustered samples'
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