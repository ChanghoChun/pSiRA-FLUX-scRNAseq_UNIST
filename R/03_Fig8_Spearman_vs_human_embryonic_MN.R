# =============================================================================
# 03_Fig8_Spearman_vs_human_embryonic_MN.R
# Figure 8: transcriptomic similarity between hiPSC-derived somatic motor
# neurons (Pulsatile, FLUX, pSiMP) and human embryonic spinal cord motor
# neurons at Carnegie stages CS12, CS14, CS17 and CS19 (Rayon et al., 2021).
# Similarity = Spearman correlation of group-averaged log-normalized
# expression, computed for all shared genes and for curated gene panels.
# =============================================================================

source("R/00_config.R")

obj <- readRDS(path_hipsc_embryo_mn_obj)

# ---- Consolidate embryonic samples into Carnegie stages ---------------------
stage_map <- c(
  "human_CS12_spinalcord_rep1" = "CS12", "human_CS12_spinalcord_rep2" = "CS12",
  "human_CS14_brachial_rep1"   = "CS14", "human_CS14_thoracic_rep1"   = "CS14",
  "human_CS17_brachial_rep1"   = "CS17", "human_CS17_brachial_rep2"   = "CS17",
  "human_CS17_thoracic_rep1"   = "CS17", "human_CS17_thoracic_rep2"   = "CS17",
  "human_CS19_brachial_rep1"   = "CS19", "human_CS19_brachial_rep2"   = "CS19",
  "human_CS19_thoracic_rep1"   = "CS19", "human_CS19_thoracic_rep2"   = "CS19"
)
obj$stage_group <- dplyr::recode(obj$orig.ident, !!!stage_map)

# Somatic MNs only (v-MNs excluded); these are the cells used in Fig. 8
obj_smn <- subset(obj, subset = new_consolidated_annotation == somatic_mn)

# ---- Average expression per group ------------------------------------------
hipsc_groups  <- c("FLUX", "pSiMP", "Pulsatile")
embryo_stages <- c("CS12", "CS14", "CS17", "CS19")

Idents(obj_smn) <- "stage_group"
expr <- GetAssayData(obj_smn, assay = "RNA", slot = "data")
avg_by_group <- lapply(levels(obj_smn), function(id) {
  rowMeans(expr[, Idents(obj_smn) == id, drop = FALSE])
}) |> setNames(levels(obj_smn))

avg_hipsc  <- avg_by_group[hipsc_groups]
avg_embryo <- avg_by_group[embryo_stages]
shared_genes <- Reduce(intersect, lapply(avg_by_group, names))

# ---- Gene panels ------------------------------------------------------------
gene_panels <- list(
  all_genes        = shared_genes,
  RA_genes         = c("RARB", "HOXA1", "CYP26A1", "STRA6", "CRABP2",
                       "ALDH1A2", "ZFP36", "DHRS3", "FABP5", "LHX1"),
  axon_guidance    = c("SEMA3A", "NTN1", "ROBO2", "ROBO3", "EFNA5",
                       "SLIT2", "SLIT3", "NRP1", "NRP2"),
  synapse_assembly = c("AGRN", "NLGN1", "BSN", "ERC1", "NRXN1", "NRXN2",
                       "SYT1", "SYT2", "STX1A", "SYN1", "SNAP25",
                       "SLC17A6", "RIMS2", "SYT9"),
  ion_channel_all  = c("SCN1A", "SCN2A", "SCN8A", "KCNA1", "KCNA2", "KCNB1",
                       "KCNQ2", "KCNQ3", "CACNA1A", "CACNA1B", "HCN1", "HCN2",
                       "CLCN3", "CLCN5", "GRIA2", "GRIN1", "GRIN2B", "CACNB4"),
  Ach_processing   = c("SYT1", "STX1A", "STX1B", "SNAP25", "VAMP2",
                       "UNC13A", "UNC13B", "STXBP1", "RIMS1", "SLC5A7",
                       "SLC18A3", "CHAT", "ACHE", "SV2A")
)

# ---- Spearman correlation ---------------------------------------------------
spearman_rho <- function(stage, group, genes) {
  genes_use <- intersect(shared_genes, genes)
  if (length(genes_use) < 2) return(NA_real_)
  cor(avg_embryo[[stage]][genes_use], avg_hipsc[[group]][genes_use],
      method = "spearman")
}

results <- expand.grid(embryo_stage = embryo_stages,
                       hipsc_group  = hipsc_groups,
                       gene_set     = names(gene_panels),
                       stringsAsFactors = FALSE) |>
  rowwise() |>
  mutate(rho = spearman_rho(embryo_stage, hipsc_group, gene_panels[[gene_set]])) |>
  ungroup()

# expected, all genes, CS12-CS19: Pulsatile 0.864 0.878 0.875 0.874;
#                                 FLUX      0.870 0.880 0.875 0.872 (Fig. 8b)
write_csv(results, file.path(results_dir, "Fig8_Spearman_hiPSC_vs_embryonic_sMN.csv"))

save_session_info("03_Fig8")
