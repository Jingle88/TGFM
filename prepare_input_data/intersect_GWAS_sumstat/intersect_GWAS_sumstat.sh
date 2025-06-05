#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J intersect_GWAS_sumstat
#BSUB -o /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/intersect_GWAS_sumstat-%J-output.log
#BSUB -e /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/intersect_GWAS_sumstat-%J-error.log 

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/2

# Set paths and variables
workdir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/scripts/intersect_GWAS_sumstat
cd $workdir

# Define input and output directories
INPUT_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/GWAS_summary_statistics"
OUTPUT_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/GWAS_summary_statistics"

# Input files
FILE1="${INPUT_DIR}/CrohnsDisease_DeLange_NatGen2017_formatted.txt.gz"
FILE2="${INPUT_DIR}/InflammatoryBowelDisease_DeLange_NatGen2017_formatted.txt.gz"
FILE3="${INPUT_DIR}/UlcerativeColitis_DeLange_NatGen2017_formatted.txt.gz"

# Output files
OUTPUT1="${OUTPUT_DIR}/intersected_CD_GWAS_sumstat.txt.gz"
OUTPUT2="${OUTPUT_DIR}/intersected_IBD_GWAS_sumstat.txt.gz"
OUTPUT3="${OUTPUT_DIR}/intersected_UC_GWAS_sumstat.txt.gz"

# Extract SNP column
zcat $FILE1 | awk 'NR>1 {print $2}' | sort > snps1.txt
zcat $FILE2 | awk 'NR>1 {print $2}' | sort > snps2.txt
zcat $FILE3 | awk 'NR>1 {print $2}' | sort > snps3.txt

# Find common SNPs
comm -12 snps1.txt snps2.txt | comm -12 - snps3.txt > common_snps.txt

# Filter original files based on common SNPs
ezcat $FILE1 | head -1 > header.txt  # Extract header
grep -w -F -f common_snps.txt <(zcat $FILE1) | cat header.txt - | gzip > $OUTPUT1
grep -w -F -f common_snps.txt <(zcat $FILE2) | cat header.txt - | gzip > $OUTPUT2
grep -w -F -f common_snps.txt <(zcat $FILE3) | cat header.txt - | gzip > $OUTPUT3

# Cleanup
rm snps1.txt snps2.txt snps3.txt common_snps.txt header.txt

echo "Intersected files created: $OUTPUT1, $OUTPUT2, $OUTPUT3"
