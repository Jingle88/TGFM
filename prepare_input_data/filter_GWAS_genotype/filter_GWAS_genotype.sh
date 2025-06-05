#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J filter_GWAS_genotype
#BSUB -o /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/filter_GWAS_genotype-%J-output.log
#BSUB -e /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/filter_GWAS_genotype-%J-error.log 

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/2

# Set paths and variables
workdir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/scripts/filter_GWAS_genotype
cd $workdir

# Define input and output directories
INPUT_GWAS_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/GWAS_summary_statistics"
INPUT_PLINK_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/genotype_data/gwas" 
OUTPUT_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/genotype_data/gwas" 

# Define input file names
PLINK_PREFIX="humancoreexome_allchr_subset_included_in_ibd_analysis"

# Define output file names
INTERSECT_VARIANT_LIST="${OUTPUT_DIR}/intersect_variants_list.txt"
FILTERED_PLINK_PREFIX="${OUTPUT_DIR}/humancoreexome_filtered_plink"

# Step 1: Extract variant IDs from GWAS summary statistics
echo "Extracting variant IDs from GWAS summary statistics..."
zcat ${INPUT_GWAS_DIR}/intersected_IBD_GWAS_sumstat.txt.gz | awk '{print $2}' | sort -V | uniq > ${INTERSECT_VARIANT_LIST}
echo "Variant list saved to ${INTERSECT_VARIANT_LIST}"

# Step 2: Extract matching variants from PLINK files
echo "Filtering PLINK files to retain only GWAS variants..."
plink2 --bfile ${INPUT_PLINK_DIR}/${PLINK_PREFIX} \
      --extract ${INTERSECT_VARIANT_LIST} \
      --make-bed \
      --out ${FILTERED_PLINK_PREFIX}

echo "Filtered PLINK files saved to ${OUTPUT_DIR}"

echo "Process completed successfully!"
