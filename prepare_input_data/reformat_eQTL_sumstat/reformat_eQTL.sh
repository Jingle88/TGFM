#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 20000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>20000] rusage[mem=20000] span[hosts=1]"
#BSUB -J reformat_eQTL_per_cellsubtype[22]
#BSUB -o /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/reformat_eQTL_per_cellsubtype-%J-output.log
#BSUB -e /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/reformat_eQTL_per_cellsubtype-%J-error.log

# Create Softpack environment
module load HGI/softpack/users/eh19/test-tgfm/1

# Define input/output directories
workdir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/scripts/reformat_eQTL
cd $workdir
input_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/IBDverse_data/eqtl_input"
output_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/eQTL_summary_statistics"

# List of chromosomes to loop through
CHR=${LSB_JOBINDEX} 

# Loop through each folder dynamically
for FOLDER in ${input_dir}/dMean__*; do
    if [ -d "$FOLDER" ]; then
        CELL_TYPE=$(basename "$FOLDER" | tr -s '_' | cut -d'_' -f2)
        
        # Check if CELL_SUBTYPE exists
        if [[ $(basename "$FOLDER") =~ dMean__${CELL_TYPE}_(.*)_ct_all ]]; then
            CELL_SUBTYPE="${BASH_REMATCH[1]}"  # This captures the part between the underscores
        else
            CELL_SUBTYPE=""  # No CELL_SUBTYPE found
        fi

        COVARIATES_FILE="${FOLDER}/OPTIM_pcs/base_output/base/Covariates.tsv"

                # Create output directory for each cell subtype
        if [ -z "$CELL_SUBTYPE" ]; then
            OUTPUT_FOLDER="${output_dir}/dMean__${CELL_TYPE}_ct_all"
        else
            OUTPUT_FOLDER="${output_dir}/dMean__${CELL_TYPE}_${CELL_SUBTYPE}_ct_all"
        fi
        mkdir -p "$OUTPUT_FOLDER"
        
        INPUT_FILE="${FOLDER}/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.${CHR}.tsv"
        OUTPUT_FILE="${OUTPUT_FOLDER}/reformatted_cis_nominal1_eqtl.${CHR}.tsv"
            
        if [ -f "$INPUT_FILE" ]; then
            python reformat_eQTL.py "$INPUT_FILE" "$COVARIATES_FILE" "$OUTPUT_FILE"
            echo "Processed: $INPUT_FILE -> $OUTPUT_FILE"
        else
            echo "Skipping: $INPUT_FILE or $COVARIATES_FILE not found."
        fi
    fi
done
