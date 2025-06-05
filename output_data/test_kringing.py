import sys
import glob
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
ggplot2_pkg = importr('ggplot2')

input_dir = "/nfs/users/nfs_j/jh59/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts"
# Load the data
for file_path in glob.glob(os.path.join(input_dir, '*_z_vec.npy')):
    base_name = os.path.basename(file_path)
    identifier = base_name.split('_z_vec.npy')[0]
    ld_file = os.path.join(input_dir, f"{identifier}_full_ld.npy")
    # Load .npy files
    z_vec = np.load(file_path, allow_pickle=True)
    gene_variant_full_ld = np.load(ld_file, allow_pickle=True)

    # Convert to R format
    z_vec_r = rpy2.robjects.FloatVector(z_vec.tolist())
    ld_mat_r = rpy2.robjects.r['as.matrix'](rpy2.robjects.r['matrix'](rpy2.robjects.FloatVector(gene_variant_full_ld.flatten(order='C')), nrow=gene_variant_full_ld.shape[0], ncol=gene_variant_full_ld.shape[1]))
    n = 59957
    # Run kriging_rss to get predicted z-scores
    print("Running kringing rss for z_vec and ld_mat in window ${base_name}")
    kriging_result = susieR_pkg.kriging_rss(z=z_vec_r, R=ld_mat_r, n=n)
    plot_data = np.array(kriging_result.rx2('plot').rx2('data'))
    z_std_diff = plot_data['z_std_diff']

    # Threshold residuals to identify good-fit indices (remove red dots)
    residual_threshold = 2.5
    good_indices = np.where(z_std_diff < residual_threshold)[0]
    print(f"Number of genetic elements in good indices: {len(good_indices)}")

    # Filter z, prior, and LD matrix
    z_vec_filtered = z_vec[good_indices]
    ld_filtered = gene_variant_full_ld[np.ix_(good_indices, good_indices)]
    len(z_vec)
    len(z_vec_filtered)
    z_vec_filtered_r = rpy2.robjects.FloatVector(z_vec_filtered.tolist())
    ld_mat_filtered_r = rpy2.robjects.r['as.matrix'](rpy2.robjects.r['matrix'](rpy2.robjects.FloatVector(ld_filtered.flatten(order='C')), nrow=ld_filtered.shape[0], ncol=ld_filtered.shape[1]))
    # run the kriging_rss again for filtered z_vec and ld
    print("Running kriging rss for filtered z_vec and ld_mat in window ${base_name}")
    kriging_result_filtered = susieR_pkg.kriging_rss(z=z_vec_filtered_r, R=ld_mat_filtered_r, n=n)
    diagnostic_plot_filtered = kriging_result_filtered.rx2('plot')
    # save the plot
    rpy2.robjects.r('''
    function(p, fname) {
        ggsave(filename=fname, plot=p, width=8, height=6, dpi=300)
    }''')(diagnostic_plot_filtered, f"diagnostic_plot_filtered_${base_name}.png")
    
