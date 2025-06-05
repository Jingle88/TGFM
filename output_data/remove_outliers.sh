#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J remove_outliers[18]
#BSUB -o /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/remove_outliers-%I-%J-output.log
#BSUB -e /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/remove_outliers-%I-%J-error.log 

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/2

# Set paths and variables
workdir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts
cd $workdir
TGFM_CODE_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts"
TGFM_INPUT_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/input_data_for_TGFM_per_chunk"
OUTPUT_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/npy_variant_info_remove_outliers"
WINDOW_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file"


# Set the jobidentifier as chunk number (loop through each chunk)
CHUNK_NUM=${LSB_JOBINDEX}
echo "Processing chunk $CHUNK_NUM..."
CHROM={1-22}

# Run the TGFM without sampling script
python ${TGFM_CODE_DIR}/generate_ld_info_remove_outliers.py \
    --trait-name "IBD"\
    --tgfm-input-data ${TGFM_INPUT_DIR}/corrected_flip_tgfm_chunk${CHUNK_NUM}_input_data_summary.txt \
    --parallel-job-identifier ${CHUNK_NUM} \
    --n-components 10 \
    --p-value-threshold 0.05 \
    --gene-tissue-pip-threshold 0.2 \
    --window-file ${WINDOW_DIR}/window_file_combined_corrected_flip.tsv \
    --output-dir ${OUTPUT_DIR}

echo "Finished processing chunk $CHUNK_NUM"
