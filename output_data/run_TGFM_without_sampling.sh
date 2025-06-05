#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J run_TGFM_without_sampling[51]
#BSUB -o /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/run_TGFM_without_sampling-%I-%J-output.log
#BSUB -e /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/run_TGFM_without_sampling-%I-%J-error.log 

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/2

# Set paths and variables
workdir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts
cd $workdir
TGFM_CODE_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/TGFM_scripts"
TGFM_INPUT_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/input_data_for_TGFM_per_chunk"
NO_SAMPLE_OUTPUT_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/tgfm_without_sampling_output"

# Set the jobidentifier as chunk number (loop through each chunk)
CHUNK_NUM=${LSB_JOBINDEX}
echo "Processing chunk $CHUNK_NUM..."

# Run the TGFM without sampling script
python ${TGFM_CODE_DIR}/run_tgfm_without_sampling.py \
    --trait-name "IBD"\
    --tgfm-input-data ${TGFM_INPUT_DIR}/corrected_flip_tgfm_chunk${CHUNK_NUM}_input_data_summary.txt \
    --parallel-job-identifier ${CHUNK_NUM} \
    --out ${NO_SAMPLE_OUTPUT_DIR}/corrected_flip \
    --n-components 10 \
    --p-value-threshold 0.05 \
    --gene-tissue-pip-threshold 0.2

echo "Finished processing chunk $CHUNK_NUM"
