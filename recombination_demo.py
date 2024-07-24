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

def find_lineages_and_tips(tree_file, tsv_file):
    tree = Phylo.read(tree_file, 'newick')
    tsv_data = pd.read_csv(tsv_file, sep='\t')
    tsv_data['clade'] = tsv_data['clade'].apply(lambda x: x.replace("auto.", "") if "auto." in x else x)
    
    lineages = {}
    for clade in tree.get_terminals():
        if clade.name and clade.name.startswith('L'):
            tips = tsv_data[tsv_data['clade'] == clade.name]['root_id'].tolist()
            lineages[clade.name] = tips
    
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
    tip = lineage_list.loc[lineage_list.Lineage == lineage_name, 'Tip'].values[0]
    sequence = sequences[tip]
    return sequence

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

def generate_recombinant_sequences(sequences, lineage_list):
    # get the number of recombinants 
    prob = 0.50
    
    lineage_num = lineage_list['Lineage'].unique().shape[0]
    print(f"all lineages num : {lineage_num}")
    recombinant_num = math.ceil((prob * lineage_num)/(1-prob)) 
    recombinant_num_low = math.ceil((0.12 * lineage_num)/(1-0.12))
    #print(f"recombinant_num:{recombinant_num}")
    # if recombinant_num == 0:
    #     recombinant_num = 1

    # decide which lineages are used to generate recombinants
    unique_lineages = lineage_list['Lineage'].unique()
    all_combinations = list(itertools.combinations(unique_lineages, 2))
    #print(f"all_combine:{all_combinations}")
    if recombinant_num >= len(all_combinations):
        re_lineage = all_combinations
        recombinant_num_low = len(re_lineage)
    else:
        re_lineage = random.sample(all_combinations, recombinant_num)

    # list to restore recombinants
    recombinants = []
    # get sequence list according to different combinations
    for tuple in re_lineage :
        l1_seq = group_sequences(lineage_list, tuple[0], sequences)
        l2_seq = group_sequences(lineage_list, tuple[1], sequences)
        group_name = f"{tuple[0]}X{tuple[1]}"

        parent1 = l1_seq
        parent2 = l2_seq
        events_range = (1,5)
        new_seq_str = simulate_recombination(str(parent1.seq), str(parent2.seq), random.choice(events_range))
        new_seq = Seq(new_seq_str)
        recombinant_id = f"{group_name}"
        recombinant = SeqRecord(new_seq, id=recombinant_id, description=f"recombination_group:high")
        recombinants.append(recombinant)
    
    return recombinants,recombinant_num_low

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
    parser.add_argument('--tsv_file', type=str, required=True, help="TSV file with tip and clade data")
    #parser.add_argument('--out', type=str, required=True, help="Output file name")
    args = parser.parse_args()
    # file name preparation
    basename = os.path.basename(args.tree).replace("_annoted.tre","") #tree/fasta name preparation
    fasta_name = os.path.join(args.fasta_dir,basename+"_seqfile.fa")

    # parse tree file
    lineage_list = find_lineages_and_tips(args.tree, args.tsv_file)
    if lineage_list['Lineage'].nunique() == 0 :
        open('empty.fa', 'w').close()
        open('empty.csv', 'w').close()
        print("Only one unique lineage found. Created empty files and exiting.")
        return

    sequences = read_sequences(fasta_name)

    # generate 'low' and 'high' recombinants
    high_recombinants,low_num = generate_recombinant_sequences(sequences,lineage_list)
    low_recombinants = random.sample(high_recombinants, low_num)
    
    # Write low recombinants
    low_fasta = f"{basename}_lowre_seqfile.fa"
    write_sequences([sequences['root']] + low_recombinants, low_fasta)

    # Write high recombinants
    high_fasta = f"{basename}_highre_seqfile.fa"
    write_sequences([sequences['root']] + high_recombinants, high_fasta)
    # # for recombinant in all_recombinants:
    # #     sequences[recombinant.id] = recombinant

    # output_sequences = [sequences['root']]  # Start with the 'root' sequence
    # output_sequences += all_recombinants   # Add all recombinant sequences
    # write_sequences(output_sequences, out_fasta)
    
    # Update lineage CSVs and save
    lineage_list = add_recombinant_info(high_recombinants, lineage_list)
    out_csv = basename+"_lineages.csv"
    lineage_list.to_csv(out_csv, index=False)

if __name__ == "__main__":
    main()
