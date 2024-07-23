import argparse
from Bio import SeqIO

def generate_depth_file(fasta_file, depth_value, output_file):
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            sequence = record.seq
            break  # Only take the first sequence (named "root")
    
    with open(output_file, "w") as out:
        for i, base in enumerate(sequence, start=1):
            out.write(f"NC_045512.2\t{i}\t{base}\t{depth_value}\n")

def main():
    parser = argparse.ArgumentParser(description="Generate depth file from fasta")
    parser.add_argument("-fasta", type=str, help="Input FASTA file containing the 'root' sequence")
    parser.add_argument("-depth", type=float, help="Depth value to be assigned to each position")
    parser.add_argument("-out", type=str, help="Output TSV file")
    
    args = parser.parse_args()
    
    generate_depth_file(args.fasta, args.depth, args.out)

if __name__ == "__main__":
    main()
