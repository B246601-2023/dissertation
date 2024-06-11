import os
import argparse
from Bio import SeqIO

def extract_and_filter_fasta_sequence(input_file, output_dir, seq_dir, target_sequence):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(seq_dir, exist_ok=True)
    
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = os.path.join(output_dir, f"{base_name}_root.fa")
    filtered_output_file = os.path.join(seq_dir, f"{base_name}.fa")
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile, open(filtered_output_file, 'w') as filtered_outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id == target_sequence:
                SeqIO.write(record, outfile, "fasta")
            if record.id.startswith('T') or record.id.startswith('root'):
                SeqIO.write(record, filtered_outfile, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract a specific sequence from a FASTA file and filter sequences.")
    parser.add_argument("--input", required=True, help="Input FASTA file.")
    parser.add_argument("--out_dir", required=True, help="Output directory for the target sequence.")
    parser.add_argument("--seq_dir", required=True, help="Output directory for the filtered sequences.")
    parser.add_argument("--target_sequence", required=True, help="Target sequence identifier.")

    args = parser.parse_args()

    extract_and_filter_fasta_sequence(args.input, args.out_dir, args.seq_dir, args.target_sequence)
