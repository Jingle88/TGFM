#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J correct_flip
#BSUB -o /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/correct_flip-%J-output.log
#BSUB -e /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/correct_flip-%J-error.log 

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/2

# Set paths and variables
workdir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/scripts/filter_GWAS_genotype
cd $workdir
FILTERED_PLINK_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/genotype_data/gwas"
FILTERED_PLINK_PREFIX="${FILTERED_PLINK_DIR}/humancoreexome_filtered_plink"
OUTPUT_PREFIX_EQTL="${FILTERED_PLINK_PREFIX}_corrected_flip_eqtl"
OUTPUT_PREFIX_GWAS="${FILTERED_PLINK_PREFIX}_corrected_flip_gwas"
OUTPUT_PREFIX="${FILTERED_PLINK_PREFIX}_corrected_flip"
EQTL_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/eQTL_summary_statistics"
GWAS_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/GWAS_summary_statistics"
GWAS_FILE="${GWAS_DIR}/intersected_IBD_GWAS_sumstat.txt.gz"
SNP_LIST_DIR="/nfs/users/nfs_j/jh59/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts/check_flip_output/check_flip_output"

'''
# Align alleles using reference, and let plink fix mismatches or strand flips for eQTL
plink2 \
  --bfile "$FILTERED_PLINK_PREFIX" \
  --update-alleles "${SNP_LIST_DIR}/flipped_snps_eqtl_final.txt" \
  --make-bed \
  --out "$OUTPUT_PREFIX_EQTL"

echo "Corrected flips and saved to: ${OUTPUT_PREFIX_EQTL}.bed/bim/fam"


# Align alleles using reference, and let plink fix mismatches or strand flips for GWAS
plink2 \
  --bfile "$FILTERED_PLINK_PREFIX" \
  --update-alleles "${SNP_LIST_DIR}/flipped_snps_gwas_final.txt" \
  --make-bed \
  --out "$OUTPUT_PREFIX_GWAS"

echo "Corrected flips and saved to: ${OUTPUT_PREFIX_GWAS}.bed/bim/fam"
'''


# Align alleles using reference, and let plink fix mismatches or strand flips for GWAS + eQTL
plink2 \
  --bfile "$FILTERED_PLINK_PREFIX" \
  --update-alleles "${SNP_LIST_DIR}/flipped_snps_combined.txt" \
  --make-bed \
  --out "$OUTPUT_PREFIX"

echo "Corrected flips and saved to: ${OUTPUT_PREFIX}.bed/bim/fam"
