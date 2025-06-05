import numpy as np 

susie_alpha = np.load("/nfs/users/nfs_j/jh59/sc-eqtl-ibd/analysis/jingling_analysis/output_data/eQTL_fm//B_10_ct_ENSG00000187010_gene_model_susie_alpha.npy")
susie_mu = np.load("/nfs/users/nfs_j/jh59/sc-eqtl-ibd/analysis/jingling_analysis/output_data/eQTL_fm//B_10_ct_ENSG00000189171_gene_model_susie_mu.npy")

# Check array properties
print("Shape:", susie_alpha.shape)
print("Data type:", susie_alpha.dtype)

# Save as CSV
np.savetxt("/nfs/users/nfs_j/jh59/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts/read_mu_alpha_npy/susie_alpha.csv", susie_alpha, delimiter=",", fmt="%.4f")
np.savetxt("/nfs/users/nfs_j/jh59/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts/read_mu_alpha_npy/susie_mu.csv", susie_mu, delimiter=",", fmt="%.4f")
