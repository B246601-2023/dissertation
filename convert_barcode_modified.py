import pandas as pd
import argparse
import shutil
import os
import subprocess
import sys
import re
import time

def basic_barcode(vcf,clade_path):
    os.makedirs("barcode_tmp", exist_ok=True)
    # cmd = f"freyja barcode-build --pb {pb_path} --path {mutation_path} --outdir barcode_tmp/ --noncl"
    # sys.stdout.flush()  
    # return_code = subprocess.run(cmd, shell=True, executable="/bin/bash",
    #                              stdout=subprocess.DEVNULL,
    #                              stderr=subprocess.PIPE)
    # if return_code.returncode != 0:
    #     print(f"Error in basic_barcode: {return_code.stderr.decode('utf-8')}", file=sys.stderr)

    # df = pd.read_csv(mutation_path, sep='\t',skiprows=3)
    # def parse_paths(df):
    #     df = df.set_index('clade')
    #     # Make sure to check with new tree versions, lineages could get trimmed.
    #     df = df.drop_duplicates(keep='last')
    #     df['from_tree_root'] = df['from_tree_root'].fillna('')
    #     df['from_tree_root'] = df['from_tree_root']\
    #     .apply(lambda x: x.replace(' ', '').strip('>').split('>'))
    #     return df
    # df = parse_paths(df)

    df = parse_tree_paths(vcf)
    df_barcodes = convert_to_barcodes(df)
    df_barcodes = reversion_checking(df_barcodes)
    df_barcodes = check_mutation_chain(df_barcodes)

    df_clade = pd.read_csv(clade_path, sep="\t")
    name_to_clade = dict(zip(df_clade['root_id'], df_clade['clade']))
    
    df_barcodes.reset_index(inplace=True)
    df_barcodes['index'] = df_barcodes['index'].map(name_to_clade)
    
    df_barcodes = df_barcodes.dropna(subset=['index'])
    df_barcodes.set_index("index",inplace=True)
    print(df_barcodes)
    df_barcodes.to_csv(os.path.join("barcode_tmp", 'usher_barcodes.csv'))
    df_barcodes.reset_index().to_feather('usher_barcodes.feather')

    return 0

# def create_lineage_path(pb_path,sample_path):
#     os.makedirs("re_tree", exist_ok=True)
#     command = f"usher -i {pb_path} -v {sample_path} -d re_tree/"
#     sys.stdout.flush()  # force python to flush
#     return_code = subprocess.run(command, shell=True, executable="/bin/bash",
#                                  stdout=subprocess.DEVNULL,
#                                  stderr=subprocess.PIPE)
#     if return_code.returncode != 0:
#         print(f"Error in basic_barcode: {return_code.stderr.decode('utf-8')}", file=sys.stderr)
#     return return_code

# def parse_tree_paths(df):
#     df = df.set_index('clade')
#     df = df.drop_duplicates(keep='last')
#     df['from_tree_root'] = df['from_tree_root'].fillna('')
    
#     results = {}
#     for index, row in df.iterrows():
#         mutations = []
#         last_node = index + ":"
        
#         # check if last_node is in the row
#         if last_node in row['from_tree_root']:
#             split_result = row['from_tree_root'].split(last_node)
#             mutations.append(split_result[-1].strip())
#             nodes = split_result[0].split('node_')[1:]
#         else:
#             nodes = row['from_tree_root'].split('node_')[1:]
        
#         for node in nodes:
#             mutations.append(node.split(':')[-1].strip())

#         results[index] = [mutations]
        
#     df_result = pd.DataFrame.from_dict(results, orient='index', columns=['from_tree_root'])
#     return df_result


def parse_tree_paths(vcf_path):
    df = pd.read_csv(vcf_path,sep='\t',skiprows=3)
    result_dict = {}
    columns_to_process = df.columns[9:]

    for col in columns_to_process:
        # 初始化空列表用于存放当前列的所有筛选出的ID
        filtered_ids = []
        # 遍历当前列的每一行
        for index, value in df[col].items():
            if value != 0 and value != '.':  # 忽略0和'.'值
                # 提取当前行的ID列的值，分割并获取对应索引的变异ID
                ids = df.at[index, 'ID'].split(',')
                value = int(value)
                if len(ids) >= (value):  # 确保索引有效
                    # 由于value是从1开始，而索引从0开始，所以需要value-1
                    filtered_ids.append(ids[value - 1])
        result_dict[col] = filtered_ids
    
    result_df = pd.DataFrame({'from_tree_root': pd.Series(result_dict)})
    # result_df = pd.DataFrame.from_dict(result_dict, orient='index', columns=['from_tree_root'])
    print(result_df)
    return result_df

def sortFun(x):
    # sort based on nuc position, ignoring nuc identities
    return int(x[1:(len(x)-1)])


def convert_to_barcodes(df):
    start_time = time.time()
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
    end_time = time.time()  # 结束时间
    print(f"convert_to_barcode function took {end_time - start_time} seconds to execute")
    return df_barcodes


# def convert_to_barcodes(df):
#     start_time = time.time()

#     # 使用字典构建初始的barcode表示，这可以避免在DataFrame中不断插入数据
#     barcode_dict = {}
#     for clade in df.index:
#         mutations = df.at[clade, 'from_tree_root']
#         # 确保mutations是可以迭代的，且如果是字符串需要正确处理
#         if isinstance(mutations, str):
#             mutations = mutations.split(',')  # 假设mutations是逗号分隔的字符串
#         for mutation in mutations:
#             if mutation not in barcode_dict:
#                 barcode_dict[mutation] = {}
#             if clade not in barcode_dict[mutation]:
#                 barcode_dict[mutation][clade] = 0
#             barcode_dict[mutation][clade] += 1

#     # 转换字典为DataFrame
#     df_barcodes = pd.DataFrame.from_dict(barcode_dict, orient='index').fillna(0).astype(int).T

#     # 确保多重突变被正确处理
#     # 预处理分割列
#     temp = {}
#     dropList = []
#     for col in df_barcodes.columns:
#         if ',' in col:
#             mutations = col.split(',')
#             for mutation in mutations:
#                 if mutation not in temp:
#                     temp[mutation] = df_barcodes[col]
#                 else:
#                     temp[mutation] += df_barcodes[col]
#             dropList.append(col)

#     # 一次性添加新列并移除不需要的列
#     df_barcodes.drop(columns=dropList, inplace=True)
#     if temp:
#         df_barcodes = pd.concat([df_barcodes, pd.DataFrame(temp)], axis=1).groupby(axis=1, level=0).sum()

#     end_time = time.time()
#     print(f"convert_to_barcode function took {end_time - start_time} seconds to execute")
#     return df_barcodes


def reversion_checking(df_barcodes):
    start_time = time.time()
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
    end_time = time.time()  # 结束时间
    print(f"reversion function took {end_time - start_time} seconds to execute")


    return df_barcodes


def identify_chains(df_barcodes):
    start_time = time.time()
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
    end_time = time.time()  # 结束时间
    print(f"identify function took {end_time - start_time} seconds to execute")

    return seq_muts

# def identify_chains(df_barcodes):
#     start_time=time.time()
#     mutation_dict = {col: col[0:len(col)-1] for col in df_barcodes.columns}
#     flip_sites = {col: col[-1] + col[1:len(col) - 1] for col in df_barcodes.columns}
    
#     seq_muts = []
#     observed_seqs = {}
    
#     for d in df_barcodes.columns:
#         base_site = mutation_dict[d]
#         flipped_site = flip_sites[d]

#         for candidate in df_barcodes.columns:
#             if base_site == mutation_dict[candidate] and flipped_site != flip_sites[candidate]:
#                 pair = sorted([d, candidate])
#                 key = list(pair + [base_site + candidate[-1]])  # 将元组转换为列表
                
#                 if tuple(key) not in observed_seqs:
#                     # Confirm observed mutation sequence
#                     if (df_barcodes[d] > 0).any() and (df_barcodes[candidate] > 0).any():
#                         observed_seqs[tuple(key)] = True
#                         seq_muts.append(key)
#     end_time = time.time()  # 结束时间
#     print(f"identify function took {end_time - start_time} seconds to execute")
#     return seq_muts


def check_mutation_chain(df_barcodes):
    start_time = time.time()
    # case when (non-reversion) mutation happens in site with existing mutation
    seq_muts = identify_chains(df_barcodes)
    print(seq_muts)
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
    end_time = time.time()  # 结束时间
    print(f"check_mutation function took {end_time - start_time} seconds to execute")

    return df_barcodes


# def check_mutation_chain(df_barcodes):
#     start_time = time.time()

#     seq_muts = identify_chains(df_barcodes)
#     while seq_muts:
#         mutations_to_add = {}
#         indices_to_update = []

#         for sm in seq_muts:
#             lin_seq_indices = df_barcodes[(df_barcodes[sm[0]] > 0) & (df_barcodes[sm[1]] > 0)].index
#             indices_to_update.extend(lin_seq_indices)
#             if sm[2] not in df_barcodes.columns:
#                 mutations_to_add[sm[2]] = lin_seq_indices

#             # Update existing mutation counts
#             df_barcodes.loc[lin_seq_indices, sm[2]] = 1
#             df_barcodes.loc[lin_seq_indices, [sm[0], sm[1]]] -= 1

#         # Add new mutations as new columns
#         for mutation, indices in mutations_to_add.items():
#             df_barcodes[mutation] = 0
#             df_barcodes.loc[indices, mutation] = 1

#         # Clean up unused mutations
#         unused = df_barcodes.columns[df_barcodes.sum() == 0]
#         df_barcodes.drop(columns=unused, inplace=True)

#         # Reversion checking and identify new chains
#         df_barcodes = reversion_checking(df_barcodes)
#         seq_muts = identify_chains(df_barcodes)

#     end_time = time.time()
#     print(f"check_mutation function took {end_time - start_time} seconds to execute")
#     return df_barcodes


def merge_barcode(basic_df, recombination_df):
    
    merged_df = pd.merge(basic_df, recombination_df, how='outer').fillna(0)
    merged_df.set_index(merged_df.columns[0], inplace=True)

    return merged_df


def delete_folder(folder_path):
    if os.path.exists(folder_path):
        shutil.rmtree(folder_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some barcodes.')
    parser.add_argument('--input', type=str, help='Input vcf file path')
    parser.add_argument('--re_vcf', type=str, help='Input re_vcf file')
    parser.add_argument('--path', type=str, help='Input mutation path file')
    args = parser.parse_args()


    outname = os.path.basename(args.re_vcf)
    outname = re.sub(r'_(lowre|highre)\.vcf$', '', outname)

    #get basic barcodes by freyja
    df_basic = basic_barcode(args.input, args.path)
    #create lineage_path file
    # create_lineage_path(args.input,args.vcf)
    
    #deal with the recombination parts
    # column_names = ["clade","from_tree_root"]
    # df = pd.read_csv("re_tree/mutation-paths.txt", sep='\t',names=column_names)
    df = parse_tree_paths(args.re_vcf)
    df_barcodes = convert_to_barcodes(df)
    df_barcodes = reversion_checking(df_barcodes)
    df_barcodes = check_mutation_chain(df_barcodes)
    df_barcodes.to_csv("re_barcode.csv")

    #merge two barcodes
    df_basic = pd.read_csv("barcode_tmp/usher_barcodes.csv")
    df_barcodes = pd.read_csv("re_barcode.csv")
    df_barcodes.rename(columns={'Unnamed: 0': 'index'}, inplace=True)
    print(df_barcodes)
    df_results = merge_barcode(df_basic, df_barcodes)

    # if "lowre" in args.vcf:
    #     suffix = "_low_barcode"
    # elif "highre" in args.vcf:
    #     suffix = "_high_barcode"
    
    df_results.to_csv(outname+".csv")
    df_results.reset_index().to_feather(outname+".feather")

    delete_folder("re_tree")
    delete_folder("barcode_tmp")
    os.remove("re_barcode.csv")
