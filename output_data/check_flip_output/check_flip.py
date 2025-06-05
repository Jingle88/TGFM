# use PLINK to compute AF from GWAS genotype data
#plink2 --bfile humancoreexome_allchr_subset_included_in_ibd_analysis --freq --out humancoreexome_allchr_subset_included_in_ibd_analysis_freq

import sys
import pandas as pd
import gzip
import os

# collection of all variants with allele flips in a dataframe
all_flipped_snps_eqtl = pd.DataFrame()
all_flipped_snps_gwas = pd.DataFrame()

def read_bim_file(bim_path):
    cols = ['CHR', 'SNP_genotype', 'CM', 'BP', 'A1', 'A2']
    return pd.read_csv(bim_path, sep="\t", header=None, names=cols)

def is_snp(a1, a2):
    return len(a1) == 1 and len(a2) == 1

def is_indel(a1, a2):
    return len(a1) > 1 or len(a2) > 1

def are_complementary(a1, a2, ref_a1, ref_a2):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return (
        comp.get(a1, None) == ref_a1 and
        comp.get(a2, None) == ref_a2
    ) or (
        comp.get(a1, None) == ref_a2 and
        comp.get(a2, None) == ref_a1
    )

def alleles_differ(a1_sumstat, a2_sumstat, a1_geno, a2_geno):
    if not (is_snp(a1_sumstat, a2_sumstat) and is_snp(a1_geno, a2_geno)):
        return False  # Skip indels
    # Check for exact match
    if (a1_sumstat == a1_geno) and (a2_sumstat == a2_geno):
        return False

    # Check for strand flip (complementary alleles)
    if are_complementary(a1_sumstat, a2_sumstat, a1_geno, a2_geno):
        return "Strand flip"

    # Check for major/minor flip (same alleles, swapped order)
    if {a1_sumstat, a2_sumstat} == {a1_geno, a2_geno}:
        return "Major/minor flip"

    # True SNP mismatch
    return "True SNP mismatch"

def check_overlap(df, genotype, mode):
    suffix = f"_{mode}"
    df['CHR'] = df['CHR'].astype(str)
    genotype['CHR'] = genotype['CHR'].astype(str)
    merged = pd.merge(df, genotype, on=['CHR', 'BP'], suffixes=(suffix, '_genotype'))

    merged['alleles_differ'] = merged.apply(
        lambda row: alleles_differ(row[f"A1{suffix}"], row[f"A2{suffix}"],
                                   row["A1_genotype"], row["A2_genotype"]),
        axis=1
    )
    # Filter only differing alleles
    diff = merged[merged['alleles_differ'].isin([
    "Strand flip",
    "True SNP mismatch",
    "Major/minor flip"])].copy()

    return diff

def check_eqtl(eqtl_file, bim_file, out_dir):
    global all_flipped_snps_eqtl

    # Extract cell subtype from filename or folder structure
    cell_subtype = os.path.basename(os.path.dirname(os.path.dirname(eqtl_file))).replace("dMean__", "")
    eqtl = pd.read_csv(eqtl_file, sep="\t")
    #rename columns to match genotype
    eqtl.rename(columns={"SNP_BP": "BP"}, inplace=True)
    eqtl['BP'] = eqtl['BP'].astype(int)

    genotype = read_bim_file(bim_file)
    result = check_overlap(eqtl, genotype, mode="eqtl")
    if result.empty:
        print(f"[INFO] No allele differences found for {cell_subtype}.")
        return
    
    flipped_snps = result[["SNP_genotype", "A1_genotype", "A2_genotype", "A1_eqtl", "A2_eqtl"]].drop_duplicates()
    output_file = os.path.join(out_dir, "flipped_snps_left_eqtl_temp.txt")
    # Append to temporary file
    flipped_snps.to_csv(output_file, mode='a', index=False, header=False, sep="\t")


def check_gwas(gwas_file, bim_file, out_dir):
    global all_flipped_snps_gwas

    col = ['CHR', 'SNP', 'BP', 'A1', 'A2', 'N', 'BETA', 'BETA_VAR']
    if gwas_file.endswith('.gz'):
        with gzip.open(gwas_file, 'rt') as f:
            gwas = pd.read_csv(f, sep="\t", header=None, names=col)
    gwas['BP'] = gwas['BP'].astype(int)
    genotype = read_bim_file(bim_file)
    
    result = check_overlap(gwas, genotype, mode="gwas")
    
    if result.empty:
        print(f"[INFO] No allele differences found for {gwas_file}.")
        return
    flipped_snps = result[["SNP_genotype", "A1_genotype", "A2_genotype", "A1_gwas", "A2_gwas"]].drop_duplicates()
    all_flipped_snps_gwas = pd.concat([all_flipped_snps_gwas, flipped_snps], ignore_index=True)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python check_flip.py [eqtl|gwas] <input_file> <bim_file>")
        sys.exit(1)

    mode = sys.argv[1]
    input_file = sys.argv[2]
    bim_file = sys.argv[3]
    out_dir = "check_flip_output"
    os.makedirs(out_dir, exist_ok=True)

    if mode == "eqtl":
        check_eqtl(input_file, bim_file, out_dir)

    elif mode == "gwas":
        check_gwas(input_file, bim_file, out_dir)
        if not all_flipped_snps_gwas.empty:
            unique_snps = all_flipped_snps_gwas.drop_duplicates()
            output_file = os.path.join(out_dir, "flipped_snps_left_gwas.txt")
            unique_snps.to_csv(output_file, index=False, header=False, sep="\t")
            print(f"[INFO] Final list of unique flipped SNPs (GWAS) saved to: {output_file}")

    else:
        print("Invalid mode. Use 'eqtl' or 'gwas'.")
