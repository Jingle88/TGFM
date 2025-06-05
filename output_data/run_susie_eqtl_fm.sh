#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J susie_eqtl_fm[21]
#BSUB -o /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/susie_eqtl_fm[21]-%J-output.log
#BSUB -e /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/susie_eqtl_fm[21]-%J-error.log 

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/2

# Set paths and variables
workdir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts
cd $workdir
TGFM_CODE_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/TGFM_scripts"  
PLINK_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/genotype_data/eQTL_genotype" 
EQTL_SUMSTAT_BASE="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/eQTL_summary_statistics" 
GWAS_SUMSTAT_BASE="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/GWAS_summary_statistics" 
OUTPUT_DIR="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/output_data/eQTL_fm" 


# list of folders containing eQTL data for all cell type clusters
CELL_TYPE_FOLDERS=(${EQTL_SUMSTAT_BASE}/dMean__*_ct_all)

# Loop through each chromosome
CHROM_NUM=${LSB_JOBINDEX} 
echo "Processing chromosome $CHROM_NUM..."

# Loop through each cell type folder
for FOLDER in "${CELL_TYPE_FOLDERS[@]}"; do
    if [[ -d "$FOLDER" ]]; then
        # Extract CELL_TYPE
        CELL_TYPE=$(basename "$FOLDER" | tr -s '_' | cut -d'_' -f2)

        # Extract CELL_SUBTYPE if available
        if [[ $(basename "$FOLDER") =~ dMean__${CELL_TYPE}_(.*)_ct_all ]]; then
            CELL_SUBTYPE="${BASH_REMATCH[1]}"
        else
            CELL_SUBTYPE=""
        fi

        # Define eQTL summary statistic file
        if [ -z "$CELL_SUBTYPE" ]; then
            EQTL_FOLDER="${EQTL_SUMSTAT_BASE}/dMean__${CELL_TYPE}_ct_all"
            FILE_PREFIX="${CELL_TYPE}_ct"
        else
            EQTL_FOLDER="${EQTL_SUMSTAT_BASE}/dMean__${CELL_TYPE}_${CELL_SUBTYPE}_ct_all"
            FILE_PREFIX="${CELL_TYPE}_${CELL_SUBTYPE}_ct"
        fi
    
        EQTL_SUMSTAT="${EQTL_FOLDER}/reformatted_cis_nominal1_eqtl.${CHROM_NUM}.tsv"
        GWAS_SUMSTAT="${GWAS_SUMSTAT_BASE}/intersected_IBD_GWAS_sumstat.txt.gz"
        QVAL_FILE="${EQTL_FOLDER}/Cis_eqtls_qval.tsv"

        EQTL_OUTPUT_STEM="${OUTPUT_DIR}/${FILE_PREFIX}"

        # Define PLINK genotype file stem
        PLINK_GENO_FILE_STEM="${PLINK_DIR}/plink_geno_chr${CHROM_NUM}_plink"

        # Run the fine-mapping script
        python "${TGFM_CODE_DIR}/susie_eqtl_fine_mapping_for_tgfm.py" \
            --eqtl-data-type SumStat \
            --chrom "$CHROM_NUM" \
            --genotype-stem "$PLINK_GENO_FILE_STEM" \
            --eqtl-sumstat "$EQTL_SUMSTAT" \
            --gwas-sumstat "$GWAS_SUMSTAT" \
            --filter-strand-ambiguous \
            --out "$EQTL_OUTPUT_STEM" \
            --qval-file "$QVAL_FILE" \
                
        echo "Completed: Chromosome $CHROM_NUM, CellType: $CELL_TYPE, CellSubtype: $CELL_SUBTYPE"
        
    fi
done


echo "All tasks completed!"
