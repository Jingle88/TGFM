import os
import numpy as np
import pandas as pd
import glob

'''
# define the .npy file path
input_dir = '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/eQTL_fm'
output_dir = '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/gene_pmces_csv'


# extract all the .npy files in the input directory
'''
'''
for file in glob.glob(os.path.join(input_dir, '*_gene_model_susie_pmces.npy')):
    #load the .npy file
    data = np.load(file, allow_pickle=True)
    #check if any file contain nan values
    if np.isnan(data).any():
        # get the number of NaN values
        num_nan = np.isnan(data).sum()
        print(f'{file} contains {num_nan} NaN values')
'''
'''
# Define parameters
target_chrom = "1"
start_pos = 220000001
end_pos = 223000001

# check the gene-tissue pairs containing variants in window 1:220000001-223000001
gene_tissue_pairs = set()

for txt_file in glob.glob(os.path.join(input_dir, '*_gene_variant_info.txt')):
    try:
        # Load .txt files
        df = pd.read_csv(txt_file, sep='\t', header=None, 
                        names=['chrom', 'variant_id', 'n', 'pos', 'alt', 'ref'],
                        dtype={'chrom': str, 'pos': int})
        
        # Filter variants in target region
        in_window = df[(df['chrom'] == target_chrom) & 
                       (df['pos'] >= start_pos) & 
                       (df['pos'] <= end_pos)]
        
        if not in_window.empty:
            # Extract gene-tissue identifier from filename
            base_name = os.path.basename(txt_file)
            identifier = base_name.split('_gene_variant_info.txt')[0]
            gene_tissue_pairs.add(identifier)
            print(f"{os.path.basename(txt_file)} contains {len(in_window)} variants in {target_chrom}:{start_pos}-{end_pos}")

    except Exception as e:
        print(f"Error processing {txt_file}: {str(e)}")
        continue

print(f"\nFound {len(gene_tissue_pairs)} gene-tissue pairs with variants in the target window")   
print('all file checked')
# convert these .npy files (containting variants in window 1:220000001-223000001) to .csv files
for identifier in gene_tissue_pairs:
    npy_file = os.path.join(input_dir, f"{identifier}_gene_model_susie_pmces.npy")
    
    if not os.path.exists(npy_file):
        print(f"Warning: Missing .npy file for {identifier}")
        continue
    
    try:
        # Load numpy array
        arr = np.load(npy_file, allow_pickle=True)
        
        # Convert to DataFrame based on dimensions
        if arr.ndim == 1:
            df = pd.DataFrame(arr, columns=['values'])
        elif arr.ndim == 2:
            df = pd.DataFrame(arr)
        else:
            df = pd.DataFrame(arr.reshape(-1, arr.shape[-1]))
        
        # Save as CSV
        output_file = os.path.join(output_dir, f"{identifier}_gene_model.csv")
        df.to_csv(output_file, index=False)
        print(f"Converted: {identifier} -> {output_file}")
        
    except Exception as e:
        print(f"Error converting {identifier}: {str(e)}")
        continue

print("\nConversion complete!")
'''
input_dir = "/nfs/users/nfs_j/jh59/sc-eqtl-ibd/analysis/jingling_analysis/output_data/tissue_prior_output"
output_dir = input_dir

# convert the .npy files to .csv files
data = np.load(os.path.join(input_dir, 'corrected_flip_tissue_specific_prior_per_window_abf_1:1-3000001.npy'), allow_pickle=True)
df = pd.DataFrame(data)
df.to_csv(os.path.join(output_dir, 'corrected_flip_tissue_specific_prior_per_window_abf_1:1-3000001.csv'), index=False)

