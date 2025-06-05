#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 20000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>20000] rusage[mem=20000] span[hosts=1]"
#BSUB -J check_flip[1-22]
#BSUB -o /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/check_flip-%J-output.log
#BSUB -e /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/check_flip-%J-error.log

# Create Softpack environment
module load HGI/softpack/users/eh19/test-tgfm/1

# Define input/output directories
workdir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts/check_flip_output
cd $workdir
GWAS_FILE="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/GWAS_summary_statistics/intersected_IBD_GWAS_sumstat.txt.gz"
EQTL_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/eQTL_summary_statistics"
HUMAN_CORE_FILE="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/genotype_data/gwas/humancoreexome_filtered_plink_corrected_flip_gwas.bim"

# Only run GWAS overlap once
#python check_flip.py gwas "$GWAS_FILE" "$HUMAN_CORE_FILE"


# define the job index
CHR=${LSB_JOBINDEX}

# loop through each cell cluster
for FOLDER in ${EQTL_DIR}/dMean__*; do
    if [ -d "$FOLDER" ]; then
        CELL_TYPE=$(basename "$FOLDER" | tr -s '_' | cut -d'_' -f2)
        
        # Check if CELL_SUBTYPE exists
        if [[ $(basename "$FOLDER") =~ dMean__${CELL_TYPE}_(.*)_ct_all ]]; then
            CELL_SUBTYPE="${BASH_REMATCH[1]}"  # This captures the part between the underscores
        else
            CELL_SUBTYPE=""  # No CELL_SUBTYPE found
        fi
        
        EQTL_FILE="${FOLDER}/reformatted_cis_nominal1_eqtl.${CHR}.tsv"
        
        if [ -f "$EQTL_FILE" ]; then
            echo "Processing eQTL for $FOLDER (chr$CHR)"
            python check_flip.py eqtl "$EQTL_FILE" "$HUMAN_CORE_FILE" 
        fi
    fi
done

if [ "$CHR" -eq 22 ]; then
    echo "Running final deduplication step..."
    python check_flip.py dedup dummy dummy
fi
