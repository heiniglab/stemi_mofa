###########################
#### Necessary libraries
########################


import scanpy as sc
import anndata as an
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scanorama
import os
import multiprocessing
import random
import time
import git
import sys
from datetime import date



#######################
#### Functions for the clustering
#######################

##############
## Execute PCA and Plot

def pca_and_plot(key_add, anndata_dict, anndata_dict_apply, n_comps_var, random_state_var):
    
    sc.set_figure_params(dpi=100, color_map = 'viridis_r')
    plt.rcParams['figure.figsize'] = [20, 5] 
    
    for key in anndata_dict:
        #### Reduce the dimensionality of the variable/ gene space
        sc.tl.pca(anndata_dict_apply[key], svd_solver = 'arpack', n_comps = n_comps_var, random_state = random_state_var) # 50 is default value

        ### Save the PCA result in anndata_dict
        anndata_dict[key].obsm[key_add + '_X_pca'] = anndata_dict_apply[key].obsm['X_pca']


        print('PCA distribution ' + key)
        sc.pl.pca_variance_ratio(anndata_dict_apply[key], log=True, n_pcs = 100)
        
        print('')
        print('Amount explained variance all PCA ' + key)
        print(sum(anndata_dict_apply[key].uns['pca']['variance_ratio']))
        
        
    ### check out amount of explained variance

    for key in anndata_dict:
        print('PCA variance ' + key)
        print('threshold 500: ' + str(sum(anndata_dict_apply[key].uns['pca']['variance_ratio'][:500])))
        print('threshold 100: ' + str(sum(anndata_dict_apply[key].uns['pca']['variance_ratio'][:100])))
        print('threshold 50: ' + str(sum(anndata_dict_apply[key].uns['pca']['variance_ratio'][:50])))
        print('threshold 30: ' + str(sum(anndata_dict_apply[key].uns['pca']['variance_ratio'][:30])))

    
    result = [anndata_dict, anndata_dict_apply]
        
    return(result)
    
### Explanation of parameters:
# anndata_dict - dictionary containing multiple anndata objects
# anndata_dict.cluster - a subseted form of the original dictionary containing the same cell indexes + keys
# n_comps_var = number of pca components
#   


###################
## Compute the neighborhood graph and cluster the data


def neighbors_and_cluster(key_add, anndata_dict, anndata_dict_apply, use_rep_var, random_state_var, n_neighbors_var, n_pcs_var, resolution =1):
    
    for key in anndata_dict:
    
    #### Compute neighborhood graph (connectivities and distances between the cells based on PCA space) - different variants + saving anndata_dict
        sc.pp.neighbors(anndata_dict_apply[key] , n_neighbors=n_neighbors_var, n_pcs=n_pcs_var, use_rep = use_rep_var ,key_added = key_add, random_state = random_state_var )
        anndata_dict[key].obsp[key_add + '_distances'] = anndata_dict_apply[key].obsp[key_add + '_distances']
        anndata_dict[key].obsp[key_add + '_connectivities'] = anndata_dict_apply[key].obsp[key_add + '_connectivities']
        print('Calculated neighborhood graph')

        ### Cluster the cells based on the neighborhood graph
        sc.tl.leiden(anndata_dict_apply[key], neighbors_key = key_add, key_added = key_add + '_cluster', random_state = random_state_var, resolution = resolution )
        anndata_dict[key].obs[key_add + '_cluster' ] = anndata_dict_apply[key].obs[key_add + '_cluster' ] 
        print('Clustered cells')

        ### Embed the neighborhood graph to plot on an UMAP
        sc.tl.umap(anndata_dict_apply[key], neighbors_key = key_add, random_state = random_state_var)  # later on plotting searches for 'umap' in anndata object - only one umap can be added

        anndata_dict[key].obsm['X_umap' + key_add ] =anndata_dict_apply[key].obsm['X_umap']  # store in anndata_dict

        print(key + key_add + ' variant finished')
        
        
    result = [anndata_dict, anndata_dict_apply]

    return(result)


