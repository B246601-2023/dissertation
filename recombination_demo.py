import random
import os
from Bio import Phylo
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import argparse
import itertools
import math

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

def simulate_recombination(parent1, parent2, num_event):
    recombination_points = sorted(random.sample(range(1, len(parent1)), num_event))
    
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

def generate_recombinant_sequences(sequences, rep, group, lineage_list):
    # get the number of recombinants 
    if group == "low":
        prob = 0.12
    elif group == "high":
        prob = 0.50
    
    #unique_lineages = lineage_list['lineage_name'].unique()
    lineage_num = lineage_list['Lineage'].unique().shape[0]
    #print(f"all lineages num : {lineage_num}")
    recombinant_num = math.floor((1 - prob) * (prob * lineage_num)) 
    #print(f"recombinant_num:{recombinant_num}")
    if recombinant_num == 0:
        recombinant_num = 1

    # decide which lineages are used to generate recombinants
    unique_lineages = lineage_list['Lineage'].unique()
    all_combinations = list(itertools.combinations(unique_lineages, 2))
    #print(f"all_combine:{all_combinations}")
    re_lineage = random.sample(all_combinations, recombinant_num)

    # list to restore recombinants
    recombinants = []
    # get sequence list according to different combinations
    for tuple in re_lineage :
        l1_seq = group_sequences(lineage_list, tuple[0], sequences)
        l2_seq = group_sequences(lineage_list, tuple[1], sequences)
        group_name = f"{tuple[0]}X{tuple[1]}_{group}"

    # generate repeated recombinants
        for i in range(rep):
            parent1 = random.sample(l1_seq, 1)[0]
            parent2 = random.sample(l2_seq, 1)[0]
            events_range = (1,5)
            new_seq_str = simulate_recombination(str(parent1.seq), str(parent2.seq), random.choice(events_range))
            new_seq = Seq(new_seq_str)
            recombinant_id = f"{group_name}_{i+1}"
            recombinant = SeqRecord(new_seq, id=recombinant_id, description=f"recombination_group:{group}")
            recombinants.append(recombinant)
    
    return recombinants

# add recombinants information to lineage lists
def add_recombinant_info(recombinants,lineage_list):
    for rec in recombinants:
        rec_id = rec.id
        lineage = rec.id.split("_")[0].strip()
        lineage_list.loc[len(lineage_list.index)] = [lineage,rec_id]
    return lineage_list


def main():
    parser = argparse.ArgumentParser(description="Generate recombinants from a phylogenetic tree.")
    parser.add_argument('--tree', type=str, required=True, help="Input tree file")
    parser.add_argument('--fasta_dir', type=str, required=True, help="Input fasta file dir")
    #parser.add_argument('--out', type=str, required=True, help="Output file name")
    args = parser.parse_args()
    # file name preparation
    basename = os.path.basename(args.tree).replace("_annoted.nh","") #tree/fasta name preparation
    fasta_name = os.path.join(args.fasta_dir,basename+"_seqfile.fa")
    out_fasta = basename+"_seqfile.fa"
    out_csv = basename+"_lineages.csv"

    # parse tree file
    lineage_list = find_lineages_and_tips(args.tree)
    sequences = read_sequences(fasta_name)

    # generate 'low' and 'high' recombinants
    low_recombinants = generate_recombinant_sequences(sequences,1,"low",lineage_list)
    high_recombinants = generate_recombinant_sequences(sequences,1,"high",lineage_list)
    all_recombinants = low_recombinants + high_recombinants
    
    # for recombinant in all_recombinants:
    #     sequences[recombinant.id] = recombinant

    output_sequences = [sequences['root']]  # Start with the 'root' sequence
    output_sequences += all_recombinants   # Add all recombinant sequences
    write_sequences(output_sequences, out_fasta)
    
    # update lineages csv and save
    lineage_list = add_recombinant_info(all_recombinants,lineage_list)
    lineage_list.to_csv(out_csv, index=False)

if __name__ == "__main__":
    main()
