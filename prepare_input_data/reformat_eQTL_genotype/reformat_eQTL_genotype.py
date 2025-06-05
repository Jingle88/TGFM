import pandas as pd
import sys


def reformat_eqtl_genotype(input_file, output_file):
    # Load eQTL genotype data
    col_names = ['chr', 'snp_id', 'pos_in_morgans', 'pos', 'eff', 'ref',]
    df = pd.read_csv(input_file, sep='\t', names=col_names)

    # Extract chromosome number from the first colum (chrx)
    df['chr'] = df['chr'].str.extract(r'(\d+)').astype(int)

    # Save as output file in output directory
    df.to_csv(output_file, sep='\t', index=False, header=False)


if __name__ == "__main__":
    if len(sys.argv) !=3:
        print("Usage: python reformat_eQTL_genotype.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    reformat_eqtl_genotype(input_file, output_file)
    print(f"Finished reformatting eQTL genotype data {input_file}")
