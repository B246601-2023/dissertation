import random
import os
from Bio import Phylo
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import argparse

# parse tree file first, to get all lineage-tips relationship

def find_lineages_and_tips(file_path):
    tree = Phylo.read(file_path, 'newick')
    lineages = {}
    for clade in tree.find_clades():
        if clade.name and clade.name.startswith('L.'):
            lineages[clade.name] = get_tips(clade)
    
    # change data to dataframe
    data = []
    for lineage, tips in lineages.items():
        for tip in tips:
            data.append([lineage, tip])
    df = pd.DataFrame(data, columns=['Lineage', 'Tip'])
    df = df.sort_values(by="Lineage").reset_index(drop=True) 
    return df

def get_tips(clade):
    tips = []
    for tip in clade.get_terminals():
        tips.append(tip.name)
    return tips

##

def read_sequences(file_name):
    seq_dict = SeqIO.to_dict(SeqIO.parse(file_name, "fasta"))
    return seq_dict

def write_sequences(sequences, file_path):
    SeqIO.write(sequences, file_path, "fasta")
##

# group sequences from one lineage, to simplize the comparision, only use "L.x" level
def group_sequences(lineage_list, lineage_name, sequences):
    tips = lineage_list.loc[lineage_list.Lineage == lineage_name, ["Tip"]]
    seq_list = []
    for value in tips["Tip"].values:
        seq_list.append(sequences[value])
    return seq_list

def simulate_recombination(parent1, parent2, recombination_rate):
    length = len(parent1)
    num_recombination_events = int(length * recombination_rate)
    
    recombination_points = sorted(random.sample(range(1, length), num_recombination_events))
    
    new_seq = ""
    switch = False
    last_point = 0
    
    for point in recombination_points:
        if switch:
            new_seq += parent2[last_point:point]
        else:
            new_seq += parent1[last_point:point]
        switch = not switch
        last_point = point
    
    new_seq += parent2[last_point:] if switch else parent1[last_point:]
    return new_seq

def generate_recombinant_sequences(sequences, recombination_rate, num_recombinants, group, lineage_list):
    recombinants = []
    l1_seq = group_sequences(lineage_list, "L.1", sequences)

    # decide l2 list according to group
    if group == "close":
        l2_seq = group_sequences(lineage_list, "L.2", sequences)
        group_name = "L.1XL.2"
    elif group in ["high", "low", "medium"]:
        Lname = lineage_list.Lineage.iloc[-1]
        l2_seq = group_sequences(lineage_list, Lname, sequences)
        group_name = f"L.1X{Lname}_{group}"

    # generate repeated recombinants
    for i in range(num_recombinants):
        parent1 = random.sample(l1_seq, 1)[0]
        parent2 = random.sample(l2_seq, 1)[0]
        new_seq_str = simulate_recombination(str(parent1.seq), str(parent2.seq), recombination_rate)
        new_seq = Seq(new_seq_str)
        recombinant_id = f"RL_{group_name}_{i+1}"
        recombinant = SeqRecord(new_seq, id=recombinant_id, description=f"recombination rate : {recombination_rate}\trecombination_group:{group}")
        recombinants.append(recombinant)
    return recombinants

# add recombinants information to lineage lists
def add_recombinant_info(recombinants,lineage_list):
    for rec in recombinants:
        rec_id = rec.id
        lineage = rec.id.split("_")[1].strip()
        lineage_list.loc[len(lineage_list.index)] = [lineage,rec_id]
    return lineage_list


def main():
    parser = argparse.ArgumentParser(description="Generate recombinants from a phylogenetic tree.")
    parser.add_argument('--tree', type=str, required=True, help="Input tree file")
    parser.add_argument('--fasta_dir', type=str, required=True, help="Input fasta file dir")
    #parser.add_argument('--out', type=str, required=True, help="Output file name")
    args = parser.parse_args()

    # 所有树的RL_close采用相同的重组率，RL_far在不同的树中使用不同的重组率，在划定区间内根据二项分布概率选择概率，并判断属于低/高哪一个区间
    recombination_rate_close = 2.0e-5
    num_recombinants = 6
    basename = os.path.basename(args.tree).replace("_annoted.nh","") #tree/fasta name preparation
    fasta_name = os.path.join(args.fasta_dir,basename+"_seqfile.fa")
    out_fasta = basename+"_seqfile.fa"
    out_csv = basename+"_lineages.csv"

    # parse tree file
    lineage_list = find_lineages_and_tips(args.tree)
    sequences = read_sequences(fasta_name)
    
    # generate 'close' recombinants
    close_recombinants = generate_recombinant_sequences(sequences, recombination_rate_close, num_recombinants, "close", lineage_list)
    
    # decide the recombination rates for 'far' group
    recombination_rates_far = [random.uniform(5.0e-7, 5.0e-5) for _ in range(3)]
    far_recombinants = []
    
    low_range = (5.0e-7, 1e-6)
    mid_range = (1e-6, 1e-5)
    high_range = (1e-5, 5.0e-4)  

    # 从每个范围内均匀分布随机抽取一个重组率
    low_rate = random.uniform(*low_range)
    mid_rate = random.uniform(*mid_range)
    high_rate = random.uniform(*high_range)
    for rate in low_rate,mid_rate,high_rate:
        if rate >= 1e-5:
            group = "high"
        elif rate <= 1e-6:
            group = "low"
        else:
            group = "medium"
        recombinants = generate_recombinant_sequences(sequences, rate, 2, group, lineage_list)
        far_recombinants.extend(recombinants)

    all_recombinants = close_recombinants + far_recombinants
    
    for recombinant in all_recombinants:
        sequences[recombinant.id] = recombinant

    write_sequences(list(sequences.values()), out_fasta)
    
    # update lineages csv and save
    lineage_list = add_recombinant_info(all_recombinants,lineage_list)
    lineage_list.to_csv(out_csv, index=False)

if __name__ == "__main__":
    main()
