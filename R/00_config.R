# =============================================================================
# 00_config.R
# Shared settings for the FLUX scRNA-seq figure analyses (Figures 5, 6 and 8).
# All scripts are run from the repository root and source this file.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)      # v5.3.0
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ---- Paths (relative to the repository root) --------------------------------
data_dir    <- "data"      # input Seurat objects (see README)
results_dir <- "results"   # all tables and plots are written here
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Seurat object of the three day-25 hiPSC-MN libraries after QC, doublet
# removal (scDblFinder), ambient RNA correction (SoupX), SCTransform
# normalization and scANVI-based cell-type label transfer.
path_hipsc_obj <- file.path(data_dir,
  "iPSC_filtered_merged_soupX_scanvilabeled_seurat_obj_juyoung.RDS")

# Seurat object combining hiPSC-derived MNs with human embryonic spinal cord
# MNs (Rayon et al., 2021; GSE171892) after label transfer and batch correction.
path_hipsc_embryo_mn_obj <- file.path(data_dir,
  "iPSC_Rayon_MN_label_transferred_batchcorrectd_seurat_obj.RDS")

# ---- Common definitions -----------------------------------------------------
# Library names in `orig.ident`:
#   "Pulsatile" = pulsatile RA                  (FASTQ Ctrl-MN; GEO RA_only)
#   "FLUX"      = sustained RA from RA-loaded pSiMPs (FASTQ RA-MN; GEO RA_at_pSiMP)
#   "pSiMP"     = pulsatile RA + blank pSiMPs   (FASTQ pSi-MN; GEO RA_plus_pSiMP)
groups_compared <- c("Pulsatile", "FLUX")   # groups shown in Figures 5 and 6

celltype_col   <- "new_consolidated_annotation"
umap_reduction <- "umap.scanvi"

mn_celltypes  <- c("Motor neurons (MN)", "Visceral motor neurons (v-MN)")
somatic_mn    <- "Motor neurons (MN)"

# Neuronal cell types used as the denominator in Fig. 5f
neuron_celltypes <- c("Motor neurons (MN)", "Visceral motor neurons (v-MN)",
                      "Pan neurons (Pan-N)", "Cholinergic neurons (ChN)",
                      "GABA/Glycinergic neurons (GABA/GlyN)",
                      "Glutamatergic neurons (GluN)", "Sensory neurons (SN)")

save_session_info <- function(tag) {
  writeLines(capture.output(sessionInfo()),
             file.path(results_dir, paste0("sessionInfo_", tag, ".txt")))
}
