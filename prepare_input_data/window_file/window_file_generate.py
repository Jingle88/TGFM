import numpy as np
import pandas as pd
import gzip
import os

# Input and output paths 
bim_path = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/data/genotypes/eQTL_genotypes_march_2024/plink_march_2025/"
ld_matrix_path = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/data/genotypes/eQTL_genotypes_march_2024/ldscores_march_2025/"
ld_output_path = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/window_file/LD_matrix_file_perwindow/"
variant_info_output_path = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/window_file/variant_info_perwindow/"
window_file_output_path = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/window_file/"

# Create output directories if they don't exist
os.makedirs(output_path, exist_ok=True)

# Initialize the .windows file content
windows_data = []

# Iterate through each chromosome (1-22)
for chr_num in range(1, 23):
    # Load .bim file to get SNP positions
    bim_file = f"{bim_path}/plink_geno_chr{chr_num}.bim"
    bim_df = pd.read_csv(bim_file, sep="\t", header=None, usecols=[0, 1, 3], names=['chr', 'snp', 'pos'])

    # Load the corresponding LD matrix
    ld_matrix_file = f"{ld_matrix_path}/ld_matrix_chr{chr_num}.ld.gz"
    with gzip.open(ld_matrix_file, "rt") as f:
        ld_matrix = np.loadtxt(f)

    # Generate 3 MB windows (overlapping by 1.5 MB)
    window_size = 3_000_000
    step_size = 1_500_000  

    start_pos = 0
    window_index = 1

    while start_pos < bim_df['pos'].max():
        end_pos = start_pos + window_size

        # Filter SNPs in this window
        window_snps = bim_df[(bim_df['pos'] >= start_pos) & (bim_df['pos'] < end_pos)]

        if not window_snps.empty:
            # Extract LD matrix slice
            idx = window_snps.index.tolist()
            ld_slice = ld_matrix[np.ix_(idx, idx)]

            # Save LD matrix and variant info
            window_name = f"{chr_num}:{start_pos}:{end_pos}"
            ld_outfile = f"{ld_output_path}/{window_name}_ld.npy"
            var_info_outfile = f"{variant_info_output_path}/{window_name}_variant_info.txt"

            np.save(ld_outfile, ld_slice)
            window_snps.to_csv(var_info_outfile, sep="\t", index=False)

            # Add to .windows file data
            windows_data.append([
                window_name,
                chr_num,
                start_pos,
                end_pos,
                ld_outfile,
                var_info_outfile
            ])

        # Move the window step forward (with overlap)
        start_pos += step_size
        window_index += 1

# Create .windows file
windows_df = pd.DataFrame(windows_data, columns=[
    "window_name", "chr", "window_start", "window_end", "ld_matrix_path", "variant_info_path"
])

windows_df.to_csv(f"{window_file_output_path}/windows_summary.txt", sep="\t", index=False)

print("Process completed successfully and window_file is generated! ")
