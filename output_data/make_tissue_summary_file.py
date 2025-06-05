import os
import re
import csv

#define the paths to input and output files
tissue_dir = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/eQTL_fm"
tissue_files = os.listdir(tissue_dir)
sample_size_dir = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/eQTL_summary_statistics"
sample_size_column_index = 7
output_file = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/tissue_summary_file.txt"


# step 1: Extract the tissue name and output stem from the gene summary files (generated in susie_eqtl_fm)
# Create a directory that hold the tissue names and their corresponding output file paths
tissue_info =  {}

for file in tissue_files:
    if not file.endswith("_gene_summary.txt"):
        continue

    match = re.match(r"(.+?)_ct_chr\d+_gene_summary.txt", file)
    if match:
        tissue_name = match.group(1)
        output_stem = os.path.join(tissue_dir, f"{tissue_name}")
        tissue_info[tissue_name] = output_stem


# step 2: extract the sample size from eQTL summmary statistics files (reformatted)
sample_size_dict = {}


for root, dirs, files in os.walk(sample_size_dir):
    match = re.search(r"dMean__([A-Za-z0-9_]+)_ct_all", root)
    if not match:
        continue

    tissue_name = match.group(1)
    for file in files:
        if re.match("reformatted_cis_nominal1_eqtl\.\d+\.tsv", file):
            full_path = os.path.join(root, file)
            try:
                with open(full_path, 'r') as f:
                    reader = csv.reader(f, delimiter='\t')
                    header = next(reader, None)
                    first_data_row = next(reader, None)
                    if first_data_row and len(first_data_row) > sample_size_column_index:
                        sample_size_value = first_data_row[sample_size_column_index].strip()
                        sample_size_dict[tissue_name] = sample_size_value
            except Exception as e:
                print(f"Error reading {full_path}: {e}")
            break  # Only need one file per tissue

def natural_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]

# step 3: write the tissue summary file
with open(output_file, 'w') as out:
    out.write("tissue_name\tsample_Size\tsusie_eqtl_output_stem\n")
    for tissue_name in sorted(tissue_info.keys(), key=natural_key):
        output_stem = tissue_info[tissue_name]
        sample_size_value = sample_size_dict.get(tissue_name, "NA")
        out.write(f"{tissue_name}\t{sample_size_value}\t{output_stem}\n")

print (f'The tissue summary file has been generated at {output_file}')
