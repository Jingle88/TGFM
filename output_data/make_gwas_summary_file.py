import os
import re
import gzip
import csv

# Define the paths to input directory and output file
trait_dir = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/GWAS_summary_statistics"
sample_size_column_index = 5
GWAS_sumstat_file = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/GWAS_summary_statistics/intersected_IBD_GWAS_sumstat.txt.gz"
output_file = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/gwas_summary_file.txt"

# Extract trait name from the GWAS summary statistics files
trait_info = {}

filename = os.path.basename(GWAS_sumstat_file)
trait_name = filename.replace("intersected_", "").replace("_GWAS_sumstat.txt.gz", "")
trait_info = {trait_name: GWAS_sumstat_file}


# Extract the sample size from GWAS summary statistics files
sample_size_dict = {}

try:
    with gzip.open(GWAS_sumstat_file, 'rt') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader, None)
        first_data_row = next(reader, None)
        if first_data_row and len(first_data_row) > sample_size_column_index:
            sample_size_value = first_data_row[sample_size_column_index].strip()
            sample_size_dict[trait_name] = sample_size_value
except Exception as e:
    print(f"Error reading {GWAS_sumstat_file}: {e}")
    sample_size_dict[trait_name] = "NA"


# Write the GWAS summary file
with open(output_file, 'w') as out:
    out.write("trait_name\tsample_size\tsummary_statistics\n")
    for trait_name in sorted(trait_info.keys()):
        file_path = trait_info[trait_name]
        sample_size_value = sample_size_dict.get(trait_name, "NA")
        out.write(f"{trait_name}\t{sample_size_value}\t{file_path}\n")

print(f"✅ GWAS summary file written to: {output_file}")
