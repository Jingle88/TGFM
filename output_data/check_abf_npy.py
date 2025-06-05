import numpy as np
abf = np.load('/nfs/users/nfs_j/jh59/sc-eqtl-ibd/analysis/jingling_analysis/output_data/tissue_prior_output/corrected_flip_tissue_specific_prior_per_window_abf_1:26000001-29000001.npy')
print("ABF shape:", abf.shape)
print("ABF NaNs?", np.any(np.isnan(abf)))
print("ABF min/max:", np.nanmin(abf), np.nanmax(abf))
