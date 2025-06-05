#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J run_tissue_prior
#BSUB -o /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/run_tissue_prior-%J-output.log
#BSUB -e /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/run_tissue_prior-%J-error.log 

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/2

# Set paths and variables
workdir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts
cd $workdir
TGFM_CODE_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/TGFM_scripts"
TiSSUE_SUMMARY_FILE="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/tissue_summary_file.txt"
JOB_FILE="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts/job_identifier_file.txt"
TGFM_WITHOUT_SAMPLING_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/tgfm_without_sampling_output"
OUTPUT_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/tissue_prior_output"
mkdir -p $OUTPUT_DIR


# Run the tissue-specific prior script
python $TGFM_CODE_DIR/run_tgfm_tissue_specific_prior.py \
    --trait-name "IBD" \
    --tissue-summary-file $TiSSUE_SUMMARY_FILE \
    --tgfm-parallel-job-identifier $JOB_FILE \
    --tgfm-without-sampling-output $TGFM_WITHOUT_SAMPLING_DIR/corrected_flip \
    --out $OUTPUT_DIR/corrected_flip \

echo "Finished processing tissue prior"
