import pandas as pd
import argparse
import shutil
import os
import subprocess
import sys

def basic_barcode(pb_path):
    os.makedirs("barcode_tmp", exist_ok=True)
    cmd = f"freyja barcode-build --pb {pb_path} --outdir barcode_tmp/ --noncl"
    sys.stdout.flush()  
    return_code = subprocess.run(cmd, shell=True, executable="/bin/bash",
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.PIPE)
    if return_code.returncode != 0:
        print(f"Error in basic_barcode: {return_code.stderr.decode('utf-8')}", file=sys.stderr)
    return return_code

def create_lineage_path(pb_path,sample_path):
    os.makedirs("re_tree", exist_ok=True)
    command = f"usher -i {pb_path} -v {sample_path} -d re_tree/"
    sys.stdout.flush()  # force python to flush
    return_code = subprocess.run(command, shell=True, executable="/bin/bash",
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.PIPE)
    if return_code.returncode != 0:
        print(f"Error in basic_barcode: {return_code.stderr.decode('utf-8')}", file=sys.stderr)
    return return_code

def parse_tree_paths(df):
    df = df.set_index('clade')
    df = df.drop_duplicates(keep='last')
    df['from_tree_root'] = df['from_tree_root'].fillna('')
    
    results = {}
    for index, row in df.iterrows():
        mutations = []
        last_node = index + ":"
        
        # check if last_node is in the row
        if last_node in row['from_tree_root']:
            split_result = row['from_tree_root'].split(last_node)
            mutations.append(split_result[-1].strip())
            nodes = split_result[0].split('node_')[1:]
        else:
            nodes = row['from_tree_root'].split('node_')[1:]
        
        for node in nodes:
            mutations.append(node.split(':')[-1].strip())

        results[index] = [mutations]
        
    df_result = pd.DataFrame.from_dict(results, orient='index', columns=['from_tree_root'])
    return df_result



def sortFun(x):
    # sort based on nuc position, ignoring nuc identities
    return int(x[1:(len(x)-1)])


def convert_to_barcodes(df):
    # builds simple barcodes, not accounting for reversions
    df_barcodes = pd.DataFrame()
    for clade in df.index:
        # sparse,binary encoding
        cladeSeries = pd.Series({c: df.loc[clade, 'from_tree_root']
                                      .count(c) for c in
                                 df.loc[clade, 'from_tree_root']}, name=clade)
        df_barcodes = pd.concat((df_barcodes, cladeSeries), axis=1)

    print('separating combined splits')
    df_barcodes = df_barcodes.T
    if '' in df_barcodes.columns:
        df_barcodes = df_barcodes.drop(columns='')
    df_barcodes = df_barcodes.fillna(0)
    temp = pd.DataFrame()
    dropList = []
    for c in df_barcodes.columns:
        # if column includes multiple mutations,
        # split into separate columns and concatenates
        if "," in c:
            for mt in c.split(","):
                if mt not in temp.columns:
                    temp = pd.concat((temp, df_barcodes[c].rename(mt)),
                                     axis=1)
                else:
                    # to handle multiple different groups with mut
                    temp[mt] += df_barcodes[c]
            dropList.append(c)
    df_barcodes = df_barcodes.drop(columns=dropList)
    df_barcodes = pd.concat((df_barcodes, temp), axis=1)
    df_barcodes = df_barcodes.groupby(axis=1, level=0).sum()
    return df_barcodes



def reversion_checking(df_barcodes):
    print('Checking for mutation pairs')
    flipPairs = [(d, d[-1] + d[1:len(d)-1] + d[0]) for d in df_barcodes.columns if (d[-1] + d[1:len(d)-1] + d[0]) in df_barcodes.columns]
    flipPairs = [list(fp) for fp in list(set(flipPairs))]

    print("Flip pairs found:", flipPairs)
    print("DataFrame shape before operation:", df_barcodes.shape)

    for fp in flipPairs:
        min_values = df_barcodes[list(fp)].min(axis=1)
        df_barcodes[fp[0]] = df_barcodes[fp[0]].subtract(min_values, axis=0)
        df_barcodes[fp[1]] = df_barcodes[fp[1]].subtract(min_values, axis=0)

    unused_columns = df_barcodes.columns[df_barcodes.sum(axis=0) == 0]
    print("Dropping unused columns:", unused_columns)
    df_barcodes = df_barcodes.drop(columns=unused_columns)

    return df_barcodes


def identify_chains(df_barcodes):

    sites = [d[0:len(d) - 1]for d in df_barcodes.columns]
    flip_sites = [d[-1] + d[1:len(d) - 1]for d in df_barcodes.columns]
    # for each mutation, find possible sequential mutations
    seq_muts = [[d, df_barcodes.columns[j], d[0:len(d) - 1] +
                 df_barcodes.columns[j][-1]]
                for i, d in enumerate(df_barcodes.columns)
                for j, d2 in enumerate(sites)
                if ((flip_sites[i] == sites[j]) and
                    (d[-1] + d[1:len(d) - 1] + d[0]) !=
                    df_barcodes.columns[j])]

    # confirm that mutation sequence is actually observed
    seq_muts = [sm for sm in seq_muts if df_barcodes[(df_barcodes[sm[0]] > 0) &
                (df_barcodes[sm[1]] > 0)].shape[0] > 0]

    mut_sites = [sortFun(sm[2]) for sm in seq_muts]
    # return only one mutation per site for each iteration
    seq_muts = [seq_muts[i] for i, ms in enumerate(mut_sites)
                if ms not in mut_sites[:i]]
    return seq_muts


def check_mutation_chain(df_barcodes):
    # case when (non-reversion) mutation happens in site with existing mutation
    seq_muts = identify_chains(df_barcodes)
    while len(seq_muts) > 0:
        # combine mutations string into single mutation
        for i, sm in enumerate(seq_muts):
            lin_seq = df_barcodes[(df_barcodes[sm[0]] > 0) &
                                  (df_barcodes[sm[1]] > 0)]
            if sm[2] not in df_barcodes.columns:
                # combination leads to new mutation
                newCol = pd.Series([(1 if dfi in lin_seq.index else 0)
                                   for dfi in df_barcodes.index], name=sm[2],
                                   index=df_barcodes.index)
                df_barcodes = pd.concat([df_barcodes, newCol], axis=1)
            else:
                # combining leads to already existing mutation
                # just add in that mutation
                df_barcodes.loc[lin_seq.index, sm[2]] = 1
            # remove constituent mutations
            df_barcodes.loc[lin_seq.index, sm[0:2]] -= 1
        # drop all unused mutations
        df_barcodes = df_barcodes.drop(
            columns=df_barcodes.columns[df_barcodes.sum(axis=0) == 0])
        # in case mutation path leads to a return to the reference.
        df_barcodes = reversion_checking(df_barcodes)
        seq_muts = identify_chains(df_barcodes)
    return df_barcodes

def merge_barcode(basic_df, recombination_df):
    merged_df = pd.merge(basic_df, recombination_df, how='outer').fillna(0)
    merged_df.rename(columns={'Unnamed: 0': 'index'}, inplace=True)
    merged_df.set_index(merged_df.columns[0], inplace=True)

    return merged_df


def delete_folder(folder_path):
    if os.path.exists(folder_path):
        shutil.rmtree(folder_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some barcodes.')
    parser.add_argument('--input', type=str, help='Input pb file path')
    parser.add_argument('--vcf', type=str, help='Input vcf file')
    parser.add_argument('--outname', type=str, help='Output file name')
    args = parser.parse_args()

    #get basic barcodes by freyja
    basic_barcode(args.input)
    #create lineage_path file
    create_lineage_path(args.input,args.vcf)
    
    #deal with the recombination parts
    column_names = ["clade","from_tree_root"]
    df = pd.read_csv("re_tree/mutation-paths.txt", sep='\t',names=column_names)
    df = parse_tree_paths(df)
    df_barcodes = convert_to_barcodes(df)
    df_barcodes = reversion_checking(df_barcodes)
    df_barcodes = check_mutation_chain(df_barcodes)
    df_barcodes.to_csv("re_barcode.csv")

    #merge two barcodes
    df_basic = pd.read_csv("barcode_tmp/usher_barcodes.csv")
    df_barcodes = pd.read_csv("re_barcode.csv")
    df_results = merge_barcode(df_basic, df_barcodes)

    df_results.to_csv(args.outname+".csv")
    df_results.reset_index().to_feather(args.outname+".feather")

    delete_folder("re_tree")
    delete_folder("barcode_tmp")
    os.remove("re_barcode.csv")
