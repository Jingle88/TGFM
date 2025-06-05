import pandas as pd 
import numpy as np
import os
import sys

# Reads covariates file and counts the number of columns (sample size).
def get_sample_size(covariates_file):
    cov_df = pd.read_csv(covariates_file, sep='\t')
    return cov_df.shape[1]

def reformat_eqtl(input_file, output_file):
    # Load data
    df = pd.read_csv(input_file, sep='\t')

    # Extract snp_bp from SNP ID
    df[['CHR', 'SNP_BP', 'A2', 'A1']] = df['variant_id'].str.extract(r'chr(\d+):(\d+):(\w+):(\w+)')
    df['SNP_BP'] = pd.to_numeric(df['SNP_BP'])

    # Compute Gene Coordinate
    df['GENE_COORD'] = df['SNP_BP'] - df['start_distance']

    # dynamically assign sample size
    df['N'] = get_sample_size(covariates_file)

    # Rename and compute new columns
    df.rename(columns={
        'phenotype_id': 'GENE',
        'variant_id': 'SNP',
        'slope': 'BETA',
    }, inplace=True)

    # Convert 'se' column to numeric (handle non-numeric values as NaN)
    df['slope_se'] = pd.to_numeric(df['slope_se'], errors='coerce')

    # Compute BETA_VAR as slope_se^2
    df['BETA_VAR'] = df['slope_se']**2

    # Select and reorder columns
    output_cols = ['GENE', 'SNP', 'CHR', 'GENE_COORD', 'SNP_BP', 'A1', 'A2', 'N', 'BETA', 'BETA_VAR']
    df = df[output_cols]

    # Save to file
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    input_file = sys.argv[1]  
    covariates_file = sys.argv[2]  
    output_file = sys.argv[3] 
    reformat_eqtl(input_file, output_file)
