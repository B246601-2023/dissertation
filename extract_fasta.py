import os
import argparse
from Bio import SeqIO

def extract_fasta_sequence(input_file, output_dir, target_sequence):
    # 获取输入文件的基名，不包括路径和后缀
    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = os.path.join(output_dir, f"{base_name}_root.fa")
    
    # 使用 Biopython 读取并筛选序列
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id == target_sequence:
                SeqIO.write(record, outfile, "fasta")
                break

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract a specific sequence from a FASTA file.")
    parser.add_argument("--input", required=True, help="Input FASTA file.")
    parser.add_argument("--out_dir", required=True, help="Output directory.")
    parser.add_argument("--target_sequence", required=True, help="Target sequence identifier.")

    args = parser.parse_args()

    extract_fasta_sequence(args.input, args.out_dir, args.target_sequence)
