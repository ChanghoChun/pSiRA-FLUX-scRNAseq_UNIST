# python/

Pre-processing, QC, scANVI reference training, label transfer and integration code (Y. Bae).

Run each file from inside `python/`. Paths are relative to it: inputs go in `../data/`, outputs land in `../results/` and `../figures/`. The `.qmd` files are R (Quarto) notebooks; the `.py` files use Python 3.11 with scvi-tools 1.1.5, scanpy 1.10.4 and anndata 0.11.3.

The notebooks were run interactively, chunk by chunk, not rendered top to bottom. Chunks marked `eval=FALSE` were run once and their outputs saved to `../results/`; later chunks read those saved files.

| Order | File | What it does | Main output |
|---|---|---|---|
| 1 | `01_preprocessing_QC.qmd` | Cell Ranger matrices to Seurat; scDblFinder doublet calls; MAD-based filtering of mitochondrial and ribosomal fractions (scuttle); SoupX ambient RNA correction; SCTransform | `iPSC_filtered_merged_soupX_seurat_obj.RDS`, `iPSC_filtered_merged_soupX_SCT_seurat_obj.RDS` |
| 2 | `02_scANVI_reference_training.py` | scANVI reference model on Rayon et al. 2021 + Andersen et al. 2023 human spinal cord data (`--epochs 400`) | `results/reference_scANVI/Andersen_Rayon_sub_gpu_400_1024_0.9_trained_2*` |
| 3a | `03_celltype_annotation.qmd` (first chunks) | Exports the three day-25 libraries (Pulsatile, FLUX, pSiMP) to AnnData for label transfer | `iPSC_filtered_merged_soupX_normalized_adata_obj.h5ad` |
| 3b | `04_scANVI_label_transfer.py` | Maps the day-25 libraries onto the reference model and predicts cell types | `results/scanvi/Andersen_Rayon_sub_tranferred_ipsc_data_umap.h5ad` |
| 3c | `03_celltype_annotation.qmd` (remaining chunks) | Adds the scANVI predictions and UMAP back to the Seurat object | `iPSC_filtered_merged_soupX_scanvilabeled_seurat_obj_juyoung.RDS` |
| 4a | `05_MN_integration.qmd` (up to "integrated") | Subsets MNs and v-MNs and merges them with Rayon et al. 2021 MNs | `results/integrated_MN/iPSC_rayon_vMN_MN_merged.h5ad` |
| 4b | `06_scANVI_MN_integration.py` | scANVI batch correction of the merged MNs (`--epochs 1000 --nlayers 2 --batchsize 128 --trainsize 0.8 --dispersion gene`) | `results/integrated_MN/iPSC_juyoung_Rayon_vMN_MN_top2000_redox_pc30_epoch1000_nlayer2_128_0.8_gene_trained_umap.h5ad` |
| 4c | `05_MN_integration.qmd` (remaining chunks) | Adds the integrated latent space and UMAP to the Seurat object | `results/integrated_MN/iPSC_Rayon_MN_label_transferred_batchcorrectd_seurat_obj.RDS` |

Reference inputs expected in `../data/reference/`:

- `Andersen_Rayon_sub.h5ad`: merged Rayon (GSE171892) and Andersen human spinal cord counts with harmonized `subtype_standardized` labels and `donor` batches.
- `GSE171892_merged_postqc_SCT_seurat_obj.RDS`: Rayon et al. 2021 Seurat object after QC.
