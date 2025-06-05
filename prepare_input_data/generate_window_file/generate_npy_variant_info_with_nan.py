import os
import numpy as np
import pandas as pd
import sys

def convert_ld_to_npy_bim_to_info(ld_matrix_dir, npy_output_dir, variant_info_output_dir, window_name):
    """Convert LD matrix file in .vcor2 format (long lines with no terminators) to .npy"""
    vcor2_file = os.path.join(ld_matrix_dir, f"{window_name}_corrected_flip.unphased.vcor2")
    vars_file = os.path.join(ld_matrix_dir, f"{window_name}_corrected_flip.unphased.vcor2.vars")
    npy_file = os.path.join(npy_output_dir, f"{window_name}_ld_corrected_flip.npy")

    """Convert PLINK .bim file to a tab-delimited variant info file"""
    bim_file = os.path.join(plink_file_dir, f"{window_name}_corrected_flip.bim")
    variant_info_file = os.path.join(variant_info_output_dir, f"{window_name}_variant_info_corrected_flip.txt")

    # Read the .vars file to determine the number of variants (N)
    with open(vars_file, "r") as f:
        variant_ids = f.read().splitlines()
    n_variants = len(variant_ids)

    # Read the .vcor2 file (long lines with no line terminators)
    with open(vcor2_file, "r") as f:
        # Read all lines (ignoring line terminators) and split into whitespace-separated elements
        elements = f.read().split()  # Split by any whitespace (spaces, tabs, newlines)
        
    # Convert to float32 array
    data = np.array(elements, dtype=np.float32)

    # Determine matrix shape
    expected_triangular = n_variants * (n_variants + 1) // 2
    expected_square = n_variants ** 2

    if len(data) == expected_triangular:
        # Case 1: Lower-triangular matrix
        matrix = np.zeros((n_variants, n_variants), dtype=np.float32)
        rows, cols = np.tril_indices(n_variants)
        matrix[rows, cols] = data
        matrix = matrix + matrix.T - np.diag(matrix.diagonal())  # Symmetrize
    elif len(data) == expected_square:
        # Case 2: Square matrix
        matrix = data.reshape(n_variants, n_variants)
    else:
        print(f"Unexpected data size for {window_name}")
        return

    # Read BIM file (tab-delimited by default)
    df = pd.read_csv(bim_file, sep="\t", header=None)
    
    # print windows containting nan values
    if np.isnan(matrix).any():
        print(f"Warning: {window_name} contains NaN values, processing removal of all rows/columns that contain nan values.")
        # Identify NaN rows/cols
        nan_rows = np.isnan(matrix).all(axis=1)
        nan_cols = np.isnan(matrix).all(axis=0)
        temp_mask = ~(nan_rows | nan_cols) # mask for rows/cols with all NaNs
        temp_matrix = matrix[temp_mask][:, temp_mask] # matrix without rows/cols with all NaNs
        temp_matrix_nan = np.isnan(temp_matrix) # extract rows and columns with any NaN in the filtered matrix (wihtout entire NaN rows/columns)
    
        temp_nan_rows = temp_matrix_nan.any(axis=1)
        temp_nan_cols = temp_matrix_nan.any(axis=0)
        mask = ~(temp_nan_rows | temp_nan_cols) # mask for rows/cols with any NaNs
        matrix_cleaned = temp_matrix[mask][:, mask]

        total_dropped = len(matrix) - len(matrix_cleaned)
        print(f"Processing {window_name}: {total_dropped} NaN-containing rows/cols found.")
        np.save(npy_file, matrix_cleaned)

        df_temp = df[temp_mask]
        df_filtered = df_temp[mask]
        df_filtered.to_csv(variant_info_file, sep="\t", index=False, header=False)
    else:
        np.save(npy_file, matrix)
        df.to_csv(variant_info_file, sep="\t", index=False, header=False)
    




def convert_bim_to_variant_info(plink_file_dir, variant_info_output_dir, window_name, temp_mask, mask):
    """Convert PLINK .bim file to a tab-delimited variant info file"""
    bim_file = os.path.join(plink_file_dir, f"{window_name}_corrected_flip.bim")
    variant_info_file = os.path.join(variant_info_output_dir, f"{window_name}_variant_info_corrected_flip.txt")

    # Read BIM file (tab-delimited by default)
    df = pd.read_csv(bim_file, sep="\t", header=None)

    # if there is NaN values found in this window, save as filtered tab-delimited text file
    if len(mask) != len(df):
        print(f"Warning: {window_name} contains NaN values, saving filtered variant info.")
        # Save as tab-delimited text file
        df_temp = df[temp_mask]
        df_filtered = df_temp[mask]
        df_filtered.to_csv(variant_info_file, sep="\t", index=False, header=False)
    else:
        df.to_csv(variant_info_file, sep="\t", index=False, header=False)




if __name__ == "__main__":
    # Parse command-line arguments
    ld_matrix_dir = sys.argv[1]
    plink_file_dir = sys.argv[2]
    npy_output_dir = sys.argv[3]
    variant_info_output_dir = sys.argv[4]
    window_name = sys.argv[5]

    convert_ld_to_npy_bim_to_info(ld_matrix_dir, npy_output_dir, variant_info_output_dir, window_name)
