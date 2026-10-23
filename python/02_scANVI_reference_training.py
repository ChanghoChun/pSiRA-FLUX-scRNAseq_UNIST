# =============================================================================
# 02_scANVI_reference_training.py
# Train the scANVI reference model by integrating the Andersen + Rayon
# spinal-cord reference datasets (the labelled reference used for downstream
# label transfer).
#
# Author: Yeji Bae and Jenny L. Smith
# Env:    Python 3.11, scvi-tools 1.1.5, scanpy 1.10.4 (GPU)
# Input:  ../data/reference/Andersen_Rayon_sub.h5ad
# Output: trained scANVI model + integrated reference in
#         ../results/reference_scANVI/
# Run:    python 02_scANVI_reference_training.py --epochs 400
# Note:   launch from python/ (paths are relative to it).
# =============================================================================
# %%
import os
import tempfile
import argparse

import scanpy as sc
import scvi
import seaborn as sns
import torch
from rich import print
# from scib_metrics.benchmark import Benchmarker

# %% 
# Argument Parsing
parser = argparse.ArgumentParser(description="scANVI integration")
parser.add_argument("--epochs", type=int, default=50, help="Number of epochs for training the model")
args = parser.parse_args()
EPOCH = args.epochs

# RUN example 
# python 02_scANVI_reference_training.py --epochs 400

# %%
# SET UP 
scvi.settings.seed = 20240802  
print("Last run with scvi-tools version:", scvi.__version__)
sc.set_figure_params(figsize=(6, 6), frameon=False)
sc.settings.figdir = "../figures/"
os.makedirs("../results/reference_scANVI", exist_ok=True)
sns.set_theme()
torch.set_float32_matmul_precision("high")
save_dir = tempfile.TemporaryDirectory()

#%config InlineBackend.print_figure_kwargs={"facecolor": "w"}
#%config InlineBackend.figure_format="retina"

# DATA 
adata_path = os.path.join("../data/reference/Andersen_Rayon_sub.h5ad")
adata = sc.read(
    adata_path
)
adata

# %% 
# PREPROCESSING 
# It is better to use 'SCANVI' integration over 'scVI' if your data is fully annotated. Since our reference data is annotated, i chose 'SCANVI'
# adata.raw = adata  # keep full dimension safe

# verbosity to max (0: only error, 1: warning, 2: info, 3:hint, 4: for debugging)
sc.settings.verbosity = 4 
sc.pp.filter_genes(adata, min_cells = 5) # filtered out 878 genes that are detected in less than 5 cells ( gpu -node killed  -> tried with --gpus 4)

sc.experimental.pp.highly_variable_genes(
    adata,
    flavor="pearson_residuals", 
    n_top_genes=7000, # should i decrease it? 
    layer="counts",
    batch_key="donor",
    subset=True,
)
adata1 = adata

# preprocessing to visualize UMAP before integration 
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50) 
sc.tl.umap(adata)
sc.pl.umap(adata, color=['dataset', 'donor','subtype_standardized_abbr'],
           save= f"Andersen_Rayon_gpu_before_integration.png")

# save
adata.write_h5ad("../results/reference_scANVI/Andersen_Rayon_sub_preprocessed.h5ad") # Q. how the file size decreased?
# %% 
adata_path = os.path.join("../results/reference_scANVI/Andersen_Rayon_sub_preprocessed.h5ad")
adata = sc.read(
    adata_path
)
# sc.pl.umap(adata, color=['dataset', 'donor','subtype_standardized_abbr'],
#            save= f"Andersen_Rayon_all_gpu_before_integration.png")


# SCANVI INTEGRATION 
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="donor") # An NVIDIA GPU may be present on this machine, but a CUDA-enabled jaxlib is not installed. Falling back to cpu
scvi.settings.num_threads = 40 
# scvi.settings.verbosity = 1 

# MODEL 
model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, dispersion="gene", gene_likelihood="nb") # n_layers: # of hidden layers  for encoder and decoder NNs, n_latent: Dimensionality of the latent space
# SCVI model with the following parameters: 
# n_hidden: 128, n_latent: 30, n_layers: 2, dropout_rate: 0.1, dispersion: gene, 
# gene_likelihood: nb, latent_distribution: normal.
# Training status: Not Trained
# Model's adata is minified?: False

scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    adata=adata,
    labels_key="subtype_standardized",
    unlabeled_category="Unknown",
)

# TRAIN
scanvi_model.train(accelerator = 'gpu', # gotta use gpu if you want to use gpu (much faster)
                   devices = 'auto',
                   max_epochs = EPOCH,
                   batch_size = 1024*4,
                   train_size = 0.8) # default 0.9
# Q. Can i visualizae epoch vs loss graph?? 
    
# SAVE
model_dir = "../results/reference_scANVI"
scanvi_model.save(dir_path = model_dir, 
           prefix = f'Andersen_Rayon_sub_gpu_{EPOCH}_1024_0.9_trained_2',  # 'Andersen_Rayon_gpu_50_1024_0.9_trained_'
           save_anndata = True, overwrite=True)
                  
# %%
# RETRIVE the LATENT SPACE
SCANVI_LATENT_KEY = "X_scANVI"
adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)
SCANVI_MDE_KEY = "X_scANVI_MDE"
adata.obsm[SCANVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCANVI_LATENT_KEY])

sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=0.3)
sc.pl.umap(adata, color=['dataset', 'donor','subtype_standardized_abbr'], save= f"Andersen_Rayon_gpu_all_{EPOCH}.png")

adata.write_h5ad(f"../results/reference_scANVI/Andersen_Rayon_sub_gpu_{EPOCH}_1024_0.9_trained_2_umap.h5ad") 
# adata_path = os.path.join("../results/reference_scANVI/Andersen_Rayon_gpu_50_1024_0.9.h5ad")
# adata = sc.read(
#     adata_path
# )

# %% 
# test 
adata_path = os.path.join("../results/reference_scANVI/Andersen_Rayon_sub_gpu_400_1024_0.9_trained_2_umap.h5ad")
adata = sc.read(
    adata_path
)
adata

SCANVI_LATENT_KEY = "X_scANVI"
SCANVI_MDE_KEY = "X_scANVI_MDE"
# adata.obsm[SCANVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCANVI_LATENT_KEY])
# %%
sc.settings.set_figure_params(figsize=(8, 8))
sc.pl.embedding(
    adata,
    basis=SCANVI_MDE_KEY,
    color=['dataset', 'donor','subtype_standardized_abbr'],
    frameon=False,
    size = 2,
    save= f"Andersen_Rayon_gpu_400_integrated_v2.pdf"
    #ncols=1
)


# sc.pl.umap(adata, color=['dataset', 'donor','subtype_standardized_abbr'])
# sc.pl.embedding(
#     adata,
#     basis="X_umap",
#     color=['dataset', 'donor','subtype_standardized_abbr'],
#     frameon=False,
#     ncols=1,
#     save= f"Andersen_Rayon_gpu_400_integrated_v2.pdf"
# )

# %%
