import os
import argparse
from Bio import SeqIO
import dendropy

def extract_and_filter_fasta_sequence(input_file, output_dir, seq_dir, tree_dir, target_sequence, n):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(seq_dir, exist_ok=True)
    os.makedirs("trees_clean", exist_ok=True)
    
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = os.path.join(output_dir, f"{base_name}_root.fa")
    filtered_output_file = os.path.join(seq_dir, f"{base_name}.fa")
    tree_file = base_name.replace('_seqfile','.tre')
    
    unique_sequences = {}  

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile, open(filtered_output_file, 'w') as filtered_outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id == target_sequence:
                SeqIO.write(record, outfile, "fasta")
            if record.id.startswith('T') or record.id.startswith('root'):
                seq_str = str(record.seq)
                if seq_str not in unique_sequences:
                    unique_sequences[seq_str] = record
                    SeqIO.write(record, filtered_outfile, "fasta")
    
    # Check if the number of unique sequences is greater than or equal to n
    if len(unique_sequences) >= n:
        # Load the tree
        tree_path = os.path.join(tree_dir,tree_file)
        print(tree_path)
        tree = dendropy.Tree.get(path=tree_path, schema='newick')
        
        # Collect the ids of unique sequences
        unique_ids = [record.id for record in unique_sequences.values()]
        
        # Retain tips that are in unique_ids
        tree.retain_taxa_with_labels(unique_ids)
        
        # Save the modified tree
        tree.write(path=os.path.join("trees_clean",tree_file), schema='newick')
        print(f"Modified tree saved to {tree_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract a specific sequence from a FASTA file and filter sequences, removing duplicates.")
    parser.add_argument("--input", required=True, help="Input FASTA file.")
    parser.add_argument("--out_dir", required=True, help="Output directory for the target sequence.")
    parser.add_argument("--seq_dir", required=True, help="Output directory for the filtered sequences.")
    parser.add_argument("--tree_dir", required=True, help="Input tree file path.")
    parser.add_argument("--n", type=int, required=True, help="Threshold for the number of unique sequences to trigger tree modification.")

    args = parser.parse_args()

    extract_and_filter_fasta_sequence(args.input, args.out_dir, args.seq_dir, args.tree_dir, "root", args.n)
