#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J reformat_eQTL_genotype[1-21]
#BSUB -o /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/reformat_eQTL_genotype[1-21]-%J-output.log
#BSUB -e /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/reformat_eQTL_genotype[1-21]-%J-error.log 

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/2

# Loop through each chromosome
CHR=${LSB_JOBINDEX}

# Set paths and variables
workdir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/scripts/reformat_eQTL_genotype
cd $workdir
input_file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/IBDverse_data/plink_genotypes/plink_geno_chr${CHR}_plink.bim"
output_file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/genotype_data/eQTL_genotype/plink_geno_chr${CHR}_plink.bim"

python reformat_eQTL_genotype.py "$input_file" "$output_file"

