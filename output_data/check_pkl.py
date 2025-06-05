import pickle
import os
import numpy as np


def inspect_pkl_file(pkl_file):
    with open(pkl_file, "rb") as f:
        data = pickle.load(f)
    
    print(f"Object type: {type(data)}")
    
    # If the object is a dictionary, list keys
    if isinstance(data, dict):
        print("\nKeys in the dictionary:")
        for key in data.keys():
            print(f"  - {key}")
            # Optionally print values (for small data)
            # print(f"    Value: {data[key]}")
    
    # If the object has attributes (e.g., a class instance)
    elif hasattr(data, "__dict__"):
        print("\nAttributes:")
        for attr in data.__dict__.keys():
            print(f"  - {attr}")
    
    return data


#tgfm_input = inspect_pkl_file("/nfs/users/nfs_j/jh59/sc-eqtl-ibd/analysis/jingling_analysis/output_data/input_data_for_TGFM/tgfm_1:120000001-123000001_tgfm_input_data.pkl")
'''
if isinstance(tgfm_input, dict) and "variants" in tgfm_input:
    variants = tgfm_input["variants"]
    print(f"\nNumber of variants in window: {len(variants)}")


# Check for pmces data
if isinstance(tgfm_input, dict) and "gene_eqtl_pmces" in tgfm_input:
    pmces = tgfm_input["gene_eqtl_pmces"]
    print("\n=== PMCES Analysis ===")
    
    # Print all PMCES values
    print("PMCES values:")
    np.set_printoptions(threshold=10)  # Show first/last 5 elements for large arrays
    print(pmces)
    
    # Check data type first
    if isinstance(pmces, np.ndarray):
        # Handle numeric arrays
        if np.issubdtype(pmces.dtype, np.number):
            print(f"\nNaN values: {np.isnan(pmces).sum()}")
            print(f"Zero values: {(pmces == 0).sum()}")
        else:
            print("\nPMCES array contains non-numeric data (dtype: {pmces.dtype})")
    else:
        print("\nPMCES is not stored in a numpy array")
        # Convert to numpy array if possible
        try:
            pmces_array = np.array(pmces, dtype=float)
            print(f"\nConverted to numeric array:")
            print(f"NaN values: {np.isnan(pmces_array).sum()}")
            print(f"Zero values: {(pmces_array == 0).sum()}")
        except ValueError:
            print("\nCould not convert PMCES to numeric array")
'''

#tgfm_trait_input = inspect_pkl_file("/nfs/users/nfs_j/jh59/sc-eqtl-ibd/analysis/jingling_analysis/output_data/input_data_for_TGFM/tgfm_1:120000001-123000001_tgfm_input_trait_data.pkl")
'''
# For tgfm_input_data.pkl
print("Gene-tissue pairs in window:", tgfm_input["gene_tissue_pairs"])
print("Variants in window:", tgfm_input["variants"])
print("PMCES in window:", tgfm_input["gene_eqtl_pmces"])

# For tgfm_trait_data.pkl
print("GWAS Beta: ", tgfm_trait_input["gwas_beta"])
print("GWAS Beta_se: ", tgfm_trait_input["gwas_beta_se"])
print("GWAS study names:", tgfm_trait_input["gwas_study_names"])
print("GWAS beta matrix shape:", tgfm_trait_input["gwas_beta"].shape)
'''


tgfm_without_sampling = inspect_pkl_file("/nfs/users/nfs_j/jh59/sc-eqtl-ibd/analysis/jingling_analysis/output_data/tgfm_without_sampling_output/corrected_flip_1:1-3000001_tgfm_no_sample_results.pkl")

# for tgfm_without_sampling pkl
print("gene-tissue pairs in window:", tgfm_without_sampling["gene_tissue_pairs"])
print("variants in window:", tgfm_without_sampling["variants"])
print("middle_gene_indices in window:", tgfm_without_sampling["middle_gene_indices"])
print("middle_variant_indices:", tgfm_without_sampling["middle_variant_indices"])
