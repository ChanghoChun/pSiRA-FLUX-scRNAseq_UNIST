# %% 
import os
import tempfile
import argparse
import numpy as np

import scanpy as sc
import scvi
import seaborn as sns
import torch
from rich import print


# %% 
# DATA
# reference 
integrated_ref_path = os.path.join("../results/reference_scANVI/Andersen_Rayon_sub_gpu_400_1024_0.9_trained_2_umap.h5ad")
adata_ref = sc.read(
    integrated_ref_path
)
adata_ref  # TODO: why do i have OPC and Oligodendrocyte progenitors? -> lets just change OPC to Oligodendrocyte progenitors for now 

# our data (only contain 'counts' layer)
adata = os.path.join("../results/iPSC_filtered_merged_soupX_normalized_adata_obj.h5ad")
adata_query = sc.read(
    adata
)
adata_query

# CHECK highly_variable_genes
adata_ref.var_names # 7000
np.unique(np.isin(adata_ref.var_names, adata_query.var_names), return_counts=True) # only TRUE

# MODEL 
scanvi_dir = "../results/reference_scANVI"  # reference model from 02_scANVI_reference_training.py
scanvi_model = scvi.model.SCANVI.load(dir_path = scanvi_dir, 
                                      prefix = "Andersen_Rayon_sub_gpu_400_1024_0.9_trained_2")

# %% 
# SETUP 
scvi.settings.seed = 20240802  
sc.set_figure_params(figsize=(6, 6
                              ), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
save_dir = tempfile.TemporaryDirectory()

sc.settings.verbosity = 1
sc.settings.figdir = "../figures/"
os.makedirs("../results/scanvi", exist_ok=True)
scvi.settings.num_threads = 40 

# %% 
# TRAIN 
# query model
scvi.model.SCANVI.prepare_query_anndata(adata_query, scanvi_model) # tutorial: adata_ref
scanvi_query = scvi.model.SCANVI.load_query_data(adata_query,
                                                 scanvi_model,
                                                 accelerator = "cpu",
                                                 device = 'auto') # TODO: UserWarning: adata.layers[counts] does not contain unnormalized count data. Are you sure this is what you want?

# train 
scanvi_query.train(
    accelerator= "cpu",
    max_epochs=200,
    batch_size = 1024*4,
    plan_kwargs={"weight_decay": 0.0}, # ensures the latent representation of reference cells will remain exactly the same
    check_val_every_n_epoch=10,
)

SCANVI_LATENT_KEY = "X_scANVI"
SCANVI_PREDICTIONS_KEY = "predictions_scanvi"
adata_query.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation()
adata_query.obs[SCANVI_PREDICTIONS_KEY] = scanvi_query.predict()

# save 
scanvi_query.save(dir_path = scanvi_dir, 
           prefix = f'Andersen_Rayon_sub_tranferred',  
           save_anndata = True, overwrite=True)
adata_query.write_h5ad(f"../results/scanvi/Andersen_Rayon_sub_tranferred_iPSC_data.h5ad") 


# %% VISUALZATION

# adata_query.obs['predictions_scanvi'] = adata_query.obs['predictions_scanvi'].replace('OPC', 'Oligodendrocyte progenitors (OPC)')
# sc.pl.umap(adata_query, color=['dataset', 'donor','subtype_standardized_abbr'])

# retrieve batch corrected normalized counts
# scanvi_norm_mat = scanvi_model.get_normalized_expression(library_size = 1e4)

# (16836, 30)  # (133469, 30)
# SCANVI_LATENT_KEY_2 = "X_scANVI_transffered"
# adata_query.obsm[SCANVI_LATENT_KEY_2] = scanvi_model.get_latent_representation()  
# adata_query.layers["scanvi_norm_transffered"] = scanvi_norm_mat

sc.pp.neighbors(adata_query, use_rep='X_scANVI')
sc.tl.leiden(adata_query) #  To achieve the future defaults please pass: flavor="igraph" and n_iterations=2.  directed must also be False to work with igraph's implementation.
sc.tl.umap(adata_query) 
sc.pl.umap(
    adata_query,
    color=['dataset', 'donor','predictions_scanvi'],
    frameon=False,
    ncols=1,
    save = '_iPSC.png'
)
adata_query.write_h5ad(f"../results/scanvi/Andersen_Rayon_sub_tranferred_ipsc_data_umap.h5ad") 

# sc.pp.scale(adata_query, max_value=10)
# sc.tl.pca(adata_query, svd_solver='arpack')
# sc.pp.neighbors(adata_query, n_neighbors=10, n_pcs=50) 
# sc.tl.umap(adata_query)

# import seaborn as sns
# custom_palette = sns.color_palette("tab20", 20)
# sc.pl.umap(adata_query, color=['donor','predictions_scanvi', "leiden"],
#            ncols=1,
#            save= f"_iPSC2.png",
#            palette=custom_palette)

# adata_query.write_h5ad(f"../results/scanvi/Andersen_Rayon_sub_tranferred_ipsc_data_umap.h5ad") 

# %% 
# VISUALIZATION CHECK
adata_path = os.path.join("../results/scanvi/Andersen_Rayon_sub_tranferred_ipsc_data_umap.h5ad")
adata_query = sc.read(
    adata_path
)
adata_query

# import seaborn as sns
# custom_palette = sns.color_palette("hsv", len(adata_query.obs['predictions_scanvi'].unique()))
# custom_palette = sns.color_palette("tab20", 20)
# sc.pl.umap(adata_query, color='predictions_scanvi', palette=custom_palette)
# sc.pl.umap(adata_query, color='leiden', palette=custom_palette, save = "_leiden.png")