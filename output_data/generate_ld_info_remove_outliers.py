import sys
import pandas as pd
import numpy as np 
import os 
import pdb
import pickle
import scipy.stats
import rpy2
import rpy2.robjects.numpy2ri as numpy2ri
import rpy2.robjects as ro
ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
susieR_pkg = importr('susieR')
import argparse


def get_credible_set_genes(phi_comp, cs_thresh):
	ordered_genes = np.argsort(-phi_comp)
	cs_genes = []
	cs_counter = 0.0
	for gene_index in ordered_genes:
		if cs_counter < cs_thresh:
			cs_genes.append(gene_index)
		cs_counter = cs_counter + phi_comp[gene_index]
	cs_genes = np.asarray(cs_genes)
	return cs_genes


def extract_uniform_log_prior_probabilities(variant_names, gene_names):
	n_snps = len(variant_names)
	n_genes = len(gene_names)
	n_ele = n_snps + n_genes

	ele_prob = 1.0/n_ele

	snp_prior = np.ones(n_snps)*ele_prob
	gene_prior = np.ones(n_genes)*ele_prob

	return np.log(snp_prior), np.log(gene_prior)


def extract_full_gene_variant_ld(standardized_eqtl_effects, variant_ld, window_name):

	expression_covariance = np.dot(np.dot(standardized_eqtl_effects, variant_ld), np.transpose(standardized_eqtl_effects)) # covariance between genes (ge_ld)
	np.fill_diagonal(expression_covariance, 1.0) # convert covariance to correlation (standardize variance to 1)


	dd = np.diag(1.0/np.sqrt(np.diag(expression_covariance))) 
	ge_ld = np.dot(np.dot(dd, expression_covariance),dd)


	gene_variant_ld = np.dot(standardized_eqtl_effects,variant_ld) # Ngenes X n_variants

	top = np.hstack((ge_ld, gene_variant_ld))
	bottom = np.hstack((np.transpose(gene_variant_ld), variant_ld))
	full_ld = np.vstack((top,bottom))
	
	return full_ld

def save_filtered_data(window_name, good_indices, ld_filtered, variant_info_file, output_dir, num_genes):
    # Load variant info
    variant_info = pd.read_csv(variant_info_file, sep='\t', header=None, 
                             names=['chr', 'variant_id', 'cm', 'pos', 'allele1', 'allele2'])
    
    # Identify variant indices in good indices
    variant_mask = np.array(good_indices) >= num_genes
    variant_positions = np.where(variant_mask)[0]
    original_variant_indices = np.array(good_indices)[variant_mask] - num_genes
    
    # Filter and save variant info
    filtered_variants = variant_info.iloc[original_variant_indices]
    variant_output_path = os.path.join(output_dir, f"{window_name}_variant_info_remove_outliers.txt")
    filtered_variants.to_csv(variant_output_path, sep='\t', header=False, index=False)
    
    # Extract and save variant-only LD
    variant_ld = ld_filtered[variant_positions,:][:, variant_positions]
    ld_output_path = os.path.join(output_dir, f"{window_name}_ld_remove_outliers.npy")
    np.save(ld_output_path, variant_ld)

def remove_outlier_in_ld_info(tgfm_data, gene_log_prior, var_log_prior, gene_variant_full_ld, n_components, gene_tissue_pip_threshold, window_name):

	# Get z-scores of each genetic element
	variant_z = tgfm_data['gwas_beta']/tgfm_data['gwas_beta_se'] # z-scores of gwas variants (including mediated and non-mediated variants)
	variant_z = np.nan_to_num(variant_z, nan=0.0) # treat z score of nan (becuase of missing GWAS sumstat varinats) as 0
	# Get z-scores of each gene-tissue pair
	new_gene_z = np.dot(tgfm_data['gene_eqtl_pmces'], variant_z)  # No need to normalize by genetic variance of gene because eqtl_pmces are standardized to have variance 1
	# Create vector of concatenated z-scores
	z_vec = np.hstack((new_gene_z,variant_z))
	# get information about the number of gene-tissue pairs in one window
	num_genes = len(tgfm_data['gene_tissue_pairs'])

	# Create concatenated vector of prior probs
	prior_probs = np.hstack((np.exp(gene_log_prior), np.exp(var_log_prior)))

	# Add warnings
	rpy2.robjects.r['options'](warn=1)

	# Convert to R format
	z_vec_r = rpy2.robjects.FloatVector(z_vec.tolist())
	full_ld_mat_r = rpy2.robjects.r['as.matrix'](rpy2.robjects.r['matrix'](rpy2.robjects.FloatVector(gene_variant_full_ld.flatten(order='C')), nrow=gene_variant_full_ld.shape[0], ncol=gene_variant_full_ld.shape[1]))

	# Run kriging_rss to get predicted z-scores
	kriging_result = susieR_pkg.kriging_rss(z=z_vec_r, R=full_ld_mat_r, n=tgfm_data['gwas_sample_size'])
	plot_data = np.array(kriging_result.rx2('plot').rx2('data'))
	z_std_diff = plot_data['z_std_diff']

	# Threshold residuals to identify good-fit indices (i.e. remove red dots)
	residual_threshold = 2.5
	good_indices = np.where(z_std_diff < residual_threshold)[0]
	bad_indices = np.where(z_std_diff >= residual_threshold)[0]

	# Filter z_vec, prior, and LD matrix + variant information file
	z_vec_filtered = z_vec[good_indices]
	prior_probs_filtered = prior_probs[good_indices]
	ld_filtered = gene_variant_full_ld[np.ix_(good_indices, good_indices)]
	
	return good_indices, ld_filtered


###########################
# Parse command line args
###########################
parser = argparse.ArgumentParser()
parser.add_argument('--window-file', default='None', type=str,
                    help='File overlapping windows')
parser.add_argument('--tgfm-input-data', default='None', type=str,
                    help='File summarizing which windows to run TGFM on')
parser.add_argument('--n-components', default=10, type=int,
                    help='Number of TGFM components')
parser.add_argument('--trait-name', default='None', type=str,
                    help='Name of GWAS trait')
parser.add_argument('--p-value-threshold', default=1e-5, type=float,
                    help='Only run TGFM on windows with at least 1 variant with p-value less than this threshold')
parser.add_argument('--gene-tissue-pip-threshold', default=0.2, type=float,
                    help='Threshold used by TGFM to run TGFM with alternative initializations')
parser.add_argument('--parallel-job-identifier', default='None', type=str,
                    help='String corresponding to name of this parallel run')
parser.add_argument('--output-dir', default='None', type=str)
args = parser.parse_args()


# First: Parse window file to create window_name -> variant_info file mapping
window_variant_info_map = {}
with open(args.window_file) as f:
    header = f.readline()  # Skip header
    for line in f:
        data = line.strip().split('\t')
        window_name = data[0]
        variant_info_file = data[5]  # 6th column contains variant_info file path
        window_variant_info_map[window_name] = variant_info_file

# Load in file containing TGFM input data
# Each row of the file corresponds to a genomic region/window to run TGFM on
# Columns correspond to tgfm input data files
tgfm_input_data = np.loadtxt(args.tgfm_input_data,dtype=str,delimiter='\t')
tgfm_input_data = tgfm_input_data[1:,:]

# Get n_windows on this run
n_windows = tgfm_input_data.shape[0]



# Loop through windows
for window_iter in range(n_windows):

	##############################
	# Extract relevent data corresponding to this window
	###############################
	data = tgfm_input_data[window_iter, :]
	# Name of window
	window_name = data[0]

	# LD file
	ld_file = data[1]
	# TGFM input data
	tgfm_input_pkl = data[2]
	# TGFM gwas input data
	tgfm_gwas_input_pkl = data[3]

	# Get corresponding variant info from window file
	try:
		variant_info_file = window_variant_info_map[window_name]
	except KeyError:
		raise ValueError(f"Missing variant info for window: {window_name} in window file")

	# Load in tgfm input data
	g = open(tgfm_input_pkl, "rb")
	tgfm_data = pickle.load(g)
	g.close()
	# Load in tgfm trait-gwas input data
	g = open(tgfm_gwas_input_pkl, "rb")
	tgfm_gwas_data = pickle.load(g)
	g.close()

	# get index corresponding to trait
	trait_index = np.where(tgfm_gwas_data['gwas_study_names'] == args.trait_name)[0][0]
	# Extract trait info for this trait
	gwas_sample_size = tgfm_gwas_data['gwas_sample_size'][trait_index]
	gwas_beta = tgfm_gwas_data['gwas_beta'][trait_index,:]
	gwas_beta_se = tgfm_gwas_data['gwas_beta_se'][trait_index,:]

	# Extract gwas p
	gwas_z = gwas_beta/gwas_beta_se
	gwas_p = scipy.stats.norm.sf(abs(gwas_z))*2.0

	# Ignore windows with no pvalues less than some threshold
	if np.min(gwas_p) > args.p_value_threshold:
		print('window skipped because of window pvalue threshold')
		t_pip.write(window_name + '\tNA\tNA\n')
		continue

	# Add gwas_beta to tgfm_data obj
	tgfm_data['gwas_beta'] = gwas_beta
	tgfm_data['gwas_beta_se'] = gwas_beta_se
	tgfm_data['gwas_sample_size'] = gwas_sample_size
	del tgfm_gwas_data


	# Load in LD
	ld_mat = np.load(ld_file)

	print('processing window:', window_name)

	# Extract full ld between genes, variants, and gene-variants
	gene_variant_full_ld = extract_full_gene_variant_ld(tgfm_data['gene_eqtl_pmces'], ld_mat, window_name)
	ld_mat


	# Extract prior probabilities of each genetic element being fine-mapped
	# Here, we use a uniform prior
	var_log_prior, gene_log_prior = extract_uniform_log_prior_probabilities(tgfm_data['variants'], tgfm_data['gene_tissue_pairs'])

	##############################
	# Run TGFM
	###############################
	good_indices, ld_filtered = remove_outlier_in_ld_info(tgfm_data, gene_log_prior, var_log_prior, gene_variant_full_ld, args.n_components, args.gene_tissue_pip_threshold, window_name)

	output_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/npy_variant_info_remove_outliers"

	# save the data with excluded outliers
	save_filtered_data(window_name, good_indices, ld_filtered, variant_info_file, output_dir, num_genes=len(tgfm_data['gene_tissue_pairs']))
