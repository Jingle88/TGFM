#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J copy_qval_file
#BSUB -o /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/copy_qval_file-%J-output.log
#BSUB -e /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/copy_qval_file-%J-error.log 

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/2

# Set paths and variables
workdir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/scripts
cd $workdir
QVAL_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/IBDverse_data/eqtl_input"
EQTL_SUMSTAT_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/eQTL_summary_statistics" 

# Loop through each folder dynamically
for FOLDER in ${QVAL_DIR}/dMean__*; do
    if [ -d "$FOLDER" ]; then
        CELL_TYPE=$(basename "$FOLDER" | tr -s '_' | cut -d'_' -f2)
        
        # Check if CELL_SUBTYPE exists
        if [[ $(basename "$FOLDER") =~ dMean__${CELL_TYPE}_(.*)_ct_all ]]; then
            CELL_SUBTYPE="${BASH_REMATCH[1]}"  # This captures the part between the underscores
        else
            CELL_SUBTYPE=""  # No CELL_SUBTYPE found
        fi

        QVAL_FILE="${FOLDER}/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv"

        if [ -z "$CELL_SUBTYPE" ]; then
            OUTPUT_FOLDER="${EQTL_SUMSTAT_DIR}/dMean__${CELL_TYPE}_ct_all"
        else
            OUTPUT_FOLDER="${EQTL_SUMSTAT_DIR}/dMean__${CELL_TYPE}_${CELL_SUBTYPE}_ct_all"
        fi

        QVAL_OUTPUT_FILE="${OUTPUT_FOLDER}/Cis_eqtls_qval.tsv"

        if [ -f "$QVAL_FILE" ]; then
            cp "$QVAL_FILE" "$QVAL_OUTPUT_FILE"
            echo "Processed: $QVAL_FILE -> $QVAL_OUTPUT_FILE"
        else
            echo "Skipping: $QVAL_FILE not found."
        fi
    fi
done
