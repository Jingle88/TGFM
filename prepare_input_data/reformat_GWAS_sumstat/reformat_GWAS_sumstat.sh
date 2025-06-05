#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J reformat_GWAS_CD_UC_IBD
#BSUB -o /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/reformat_GWAS_CD_UC_IBD-%J-output.log
#BSUB -e /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/reformat_GWAS_CD_UC_IBD-%J-error.log

set -x  # Debug mode: prints each command before execution
set -e  # Exit on first error

echo "Job started on $(hostname) at $(date)"
module load HGI/softpack/users/eh19/test-tgfm/1
which python
python --version

input_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/core_analysis_output/IBDverse_multi-tissue_eQTL_project/IBDverse_coloc/gwas_files"
output_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/GWAS_summary_statistics"

declare -A diseases
diseases=(
    ["CrohnsDisease"]="40266"
    ["UlcerativeColitis"]="45975"
    ["InflammatoryBowelDisease"]="59957"
)

for disease in "${!diseases[@]}"; do
    echo "Processing $disease..."
    input_file="$input_dir/${disease}_DeLange_NatGen2017.txt.gz"
    output_file="$output_dir/${disease}_DeLange_NatGen2017_formatted.txt.gz"
    sample_size="${diseases[$disease]}"

    if [ ! -f "$input_file" ]; then
        echo "Error: Input file $input_file does not exist!"
        exit 1
    fi

    python reformat_GWAS.py --input_file "$input_file" --output_file "$output_file" --sample_size "$sample_size" || { echo "Python script failed"; exit 1; }
done

echo "All files processed successfully."
