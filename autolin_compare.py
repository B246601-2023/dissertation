import pandas as pd
import argparse

def handle_data(input_file, output_file):
    """
    Reads a CSV file, processes it by removing 'auto.' prefixes, and saves it to another file.
    
    Parameters:
        input_file (str): Path to the input CSV file.
        output_file (str): Path to the output CSV file where the processed data will be saved.
    """
    columns = [0, 1]
    df = pd.read_csv(input_file, sep='\t', usecols=columns)
    df['clade'] = df['clade'].apply(lambda x: x.replace("auto.", "")).sort_values
    # can apply to plots directly
    df.to_csv(output_file, index=False)

def main():
 
    parser = argparse.ArgumentParser(description='Process and save data.')
    parser.add_argument('--input', type=str, help='Input file path for auto annotations')
    parser.add_argument('--output', type=str, help='Output file path for modified auto annotations')
    args = parser.parse_args()

    handle_data(args.input, args.output)
    # Handle pango annotations (if needed, can use handle_data similarly)
    # handle_data(args.pango_input, "path_to_save_pango_output.csv")

if __name__ == "__main__":
    main()

