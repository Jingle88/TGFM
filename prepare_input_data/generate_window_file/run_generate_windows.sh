#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J generate_windows[1-22]
#BSUB -o /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/generate_windows[1-22]-%J-output.log 
#BSUB -e /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/generate_windows[1-22]-%J-error.log

# Create Softpack environment
module load HGI/softpack/users/eh19/test-tgfm/1

# Set paths and variables
workdir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/scripts/generate_window_file
cd $workdir

# Define input PLINK file directory and output directory
PLINK_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/genotype_data/gwas"
OUTPUT_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/windows"


# Apply the Python script to generate windows for each chromosome

# List of chromosomes to loop through
CHROM=${LSB_JOBINDEX} 

BIM_FILE="${PLINK_DIR}/humancoreexome_filtered_plink_corrected_flip.bim"
OUTPUT_FILE="${OUTPUT_DIR}/windows_chr${CHROM}_corrected_flip.tsv"

# Run the Python script
python generate_windows.py "$BIM_FILE" "$OUTPUT_FILE" "$CHROM"


echo "All chromosome window files generated successfully!"
