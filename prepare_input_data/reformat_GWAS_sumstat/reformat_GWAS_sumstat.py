import pandas as pd
import gzip
import argparse

def reformat_GWAS(input_file, output_file, sample_size):
    # Load data
    col_names = ['RSid', 'Chr', 'Pos', 'Eff_allele', 'MAF', 'pval', 'beta', 'OR', 'log_OR', 'se', 'z.score', 'Disease', 'PubmedID', 'used_file']
    df = pd.read_csv(input_file, sep='\t', names=col_names, skiprows=1)

    # Convert 'se' column to numeric (handle non-numeric values as NaN)
    df['se'] = pd.to_numeric(df['se'], errors='coerce')

    # Extract reference allele from variant_id (e.g., chr1:7773966:G:A -> G)
    df['A2'] = df['RSid'].apply(lambda x: x.split(':')[2] if len(x.split(':')) >= 3 else None)

    # Create output format columns
    df['CHR'] = df['Chr']
    df['SNP'] = df['RSid']
    df['BP'] = df['Pos']
    df['A1'] = df['Eff_allele'] 
    df['N'] = sample_size
    df['BETA'] = df['beta']
    df['BETA_VAR'] = df['se'] ** 2

    # Select and reorder relevant columns
    output_cols = ['CHR', 'SNP', 'BP', 'A1', 'A2', 'N', 'BETA', 'BETA_VAR']
    df_out = df[output_cols]

    # Save to compressed .txt.gz file
    with gzip.open(output_file, 'wt') as f:
        df_out.to_csv(f, sep='\t', index=False)


 
def parse_options():
    parser = argparse.ArgumentParser(
        description="reformat_GWAS_for_UC_DC_IBD"
   )
    parser.add_argument("--input_file","--input_file", required=True, help="Input .txt.gz file containing GWAS data.")
    parser.add_argument("--output_file","--output_file", required=True, help="Output .txt.gz file with reformatted data.")
    parser.add_argument("--sample_size","--sample_size", type=int, required=True, help="Sample size for the GWAS study.")
       
    return parser.parse_args()


def main():
    inherited_options = parse_options()
    input_file = inherited_options.input_file
    output_file = inherited_options.output_file
    sample_size = inherited_options.sample_size

    print(f":) running the code with input file {input_file} and output_file {output_file} and sample_size  {sample_size}")
    print("Parsed args")

    reformat_GWAS(input_file, output_file, sample_size)
    print("finished!")

# Execute
if __name__ == '__main__':
   main()
