# %% 
import os
import tempfile
import argparse

import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns
import torch
from rich import print



# %% 
# Argument Parsing (to test out the optimal hyperparameters)
parser = argparse.ArgumentParser(description="scANVI integration")
parser.add_argument("--topgenes", type=int, default=2000, help="Number of top variable genes")
parser.add_argument("--pcnums", type=int, default=30, help="Number of PC components to be used for neighbors computation")
parser.add_argument("--epochs", type=int, default=50, help="Number of epochs for training the model")
parser.add_argument("--nlayers", type=int, default=1, help="Number of hidden layers for encoder and decoder NNs")
parser.add_argument("--batchsize", type=int, default=128, help="Assign the batch size when training the model")
parser.add_argument("--trainsize", type=float, default=0.9, help="training set percentage")
parser.add_argument("--dispersion", type=str, default="gene", help="e.g gene, gene-batch")

# RUN example 
# python 06_scANVI_MN_integration.py --epochs 1000 --nlayers 2 --batchsize 128 --trainsize 0.8 --dispersion gene
args = parser.parse_args()
TOPgenes = args.topgenes
PCnum = args.pcnums
EPOCH = args.epochs
Nlayers = args.nlayers
BATCHsize = args.batchsize
TRAINsize = args.trainsize
DISPERSION = args.dispersion

# TOPgenes = 2000
# PCnum = 30
# EPOCH = 1000
# Nlayers = 2
# BATCHsize = 128
# TRAINsize = 0.8
# DISPERSION = "gene"

PREFIX = f"iPSC_juyoung_Rayon_vMN_MN_top{TOPgenes}_redox_pc{PCnum}_epoch{EPOCH}_nlayer{Nlayers}_{BATCHsize}_{TRAINsize}_{DISPERSION}_trained"
OUTPUT = "../results/integrated_MN/"
# %%
# SET UP 
scvi.settings.seed = 20240802  
print("Last run with scvi-tools version:", scvi.__version__)
sc.set_figure_params(figsize=(6, 6), frameon=False)
sc.settings.figdir = "../figures/"
sns.set_theme()
torch.set_float32_matmul_precision("high")
save_dir = tempfile.TemporaryDirectory()
# verbosity to max (0: only error, 1: warning, 2: info, 3:hint, 4: for debugging)
sc.settings.verbosity = 4 

# DATA 
adata_path = os.path.join("../results/integrated_MN/iPSC_rayon_vMN_MN_merged.h5ad")
adata = sc.read(adata_path)

# adata = adata[adata.obs["new_consolidated_annotation"].isin(["Motor neurons (MN)", "Visceral motor neurons (v-MN)"])]
# %% 
# PREPROCESSING 
# It is better to use 'SCANVI' integration over 'scVI' if your data is fully annotated. Since our reference data is annotated, i chose 'SCANVI'
sc.pp.filter_genes(adata, min_cells = 5) # filtered out 7721 genes that are detected in less than 5 cells,,, should i decrease that?

sc.experimental.pp.highly_variable_genes(
    adata,
    flavor="pearson_residuals", 
    n_top_genes=TOPgenes, # MN (previously used only 2000 genes)
    batch_key="sample",
    subset=False, # True -> False to keep all genes
)

# preprocessing to visualize UMAP before integration 
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=5, n_pcs=PCnum) 
sc.tl.umap(adata)
sc.pl.umap(adata, color=['orig.ident', 'sample', 'new_consolidated_annotation'], ncols=1)

# SCANVI INTEGRATION 
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample") 
scvi.settings.num_threads = 40 

# MODEL 
# SCVI model with the following parameters: https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.html#scvi-model-scvi
# generally work well parameters from https://docs.scvi-tools.org/en/1.0.0/tutorials/notebooks/harmonization.html#integration-with-scvi 
# n_hidden: 128, n_latent: 30, n_layers: 2, dropout_rate: 0.1, dispersion: gene, gene_likelihood: 'nb', latent_distribution: 'normal'
model = scvi.model.SCVI(adata, n_layers=Nlayers, n_latent=30, 
                        dispersion=DISPERSION, # dafulat = "gene" 
                        gene_likelihood="nb")  # defulat = "zinb" - zero-inflated NB distribution

scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    adata=adata,
    labels_key="new_consolidated_annotation",
    unlabeled_category="Unknown",
)

# TRAIN
scanvi_model.train(accelerator = 'gpu', # gotta use gpu if you want to use gpu (much faster)
                   devices = 'auto',
                   max_epochs = EPOCH,
                   batch_size = BATCHsize,       # dafault = 128 
                   train_size = TRAINsize)       # default 0.9
    
# SAVE
model_dir = OUTPUT 
if not os.path.exists(model_dir):
    os.makedirs(model_dir)
scanvi_model.save(dir_path = model_dir, 
           prefix = PREFIX, save_anndata = True, overwrite=True)
                  
# RETRIVE the LATENT SPACE
SCANVI_LATENT_KEY = "X_scANVI"
SCANVI_MDE_KEY = "X_scANVI_MDE"
adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)
adata.obsm[SCANVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCANVI_LATENT_KEY])

sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=0.3)
sc.pl.umap(adata, color=['orig.ident', 'sample','new_consolidated_annotation'], ncols = 1)
sc.pl.umap(adata, color=['orig.ident', 'sample','new_consolidated_annotation'], ncols=1, save= f"{PREFIX}.png")
adata.write_h5ad(f"{OUTPUT}{PREFIX}_umap.h5ad") 

# # %% 
# loss 
import matplotlib.pyplot as plt

# Extract training loss
train_loss = scanvi_model.history["elbo_train"]

# Get the number of epochs
epochs = range(1, len(train_loss) + 1)

# Plot ELBO training loss
plt.figure(figsize=(8, 6))
plt.plot(epochs, train_loss, label="Training Loss (ELBO)", color="blue")
plt.xlabel("Epochs")
plt.ylabel("ELBO Loss")
plt.title("Training Loss Over Epochs")
plt.legend()
plt.savefig(f"{OUTPUT}{PREFIX}_trained_ELBOtrainingloss.png", dpi=300, bbox_inches='tight')
plt.show()

# %% 
# visualization 
# PREFIX = "iPSC_rayon_vMN_MN_top2000_redox_pc30_epoch1000_nlayer2_128_0.8_gene_trained"
# adata_path = os.path.join(f"results/integrated_MN/{PREFIX}_umap.h5ad")
# adata = sc.read(
#     adata_path
# )

# SCANVI_LATENT_KEY = "X_scANVI"
# SCANVI_MDE_KEY = "X_scANVI_MDE"
# # adata.obsm[SCANVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCANVI_LATENT_KEY])
# sc.pl.embedding(
#     adata,
#     basis=SCANVI_MDE_KEY,
#     color=['orig.ident', 'sample','new_consolidated_annotation'],
#     frameon=False,
#     ncols=1,
#     save= f"{PREFIX}_umap.png"
# )
# %% 
##############
# DE testing #
##############
# model_dir = "results/integrated_MN"
# scanvi_model = scvi.model.SCANVI.load(dir_path=model_dir, prefix= f'{PREFIX}')

# adata.obs["sample"].value_counts()

# cell_idx1 = adata.obs["sample"] == "iPSC"
# cell_idx2 = adata.obs["sample"] == "Rayon"
# de_change = scanvi_model.differential_expression(idx1=cell_idx1, idx2=cell_idx2)
# de_change.to_csv(f'results/integrated_MN/DEGs/{PREFIX}_trained_DEGs.csv', index=True)
