# =============================================================================
# 02_Fig6b_DEG_FLUX_vs_Pulsatile.R
# Differentially expressed genes (DEGs) between FLUX and Pulsatile motor
# neurons (MN + v-MN clusters), Seurat FindMarkers on SCT-normalized data.
#   Fig. 6b  Volcano plot: 240 DEGs, 79 higher in FLUX and 161 higher in
#            Pulsatile (|log2FC| >= 0.5, adjusted p < 0.05)
# =============================================================================

source("R/00_config.R")

obj_all <- readRDS(path_hipsc_obj)
obj_mn <- subset(obj_all, subset = orig.ident %in% groups_compared &
                                   new_consolidated_annotation %in% mn_celltypes)
rm(obj_all)
table(obj_mn$orig.ident)

# ---- DEG analysis: FLUX vs Pulsatile ----------------------------------------
Idents(obj_mn) <- "orig.ident"
obj_mn <- PrepSCTFindMarkers(obj_mn)

deg <- FindMarkers(obj_mn, ident.1 = "FLUX", ident.2 = "Pulsatile",
                   assay = "SCT", slot = "data",
                   logfc.threshold = 0.5, min.pct = 0.25) %>%
  rownames_to_column("gene") %>%
  mutate(comparison = "FLUX vs Pulsatile", .before = 1)
# avg_log2FC > 0: higher in FLUX; avg_log2FC < 0: higher in Pulsatile
write_csv(deg, file.path(results_dir, "Fig6b_DEG_FLUX_vs_Pulsatile_MN.csv"))

deg_sig <- deg %>% filter(p_val_adj < 0.05)
count(deg_sig, direction = ifelse(avg_log2FC > 0, "Higher in FLUX",
                                                   "Higher in Pulsatile"))
# expected: Higher in FLUX 79, Higher in Pulsatile 161

# ---- Fig. 6b: Volcano plot --------------------------------------------------
volcano_df <- deg_sig %>%
  mutate(neglog10p = -log10(p_val_adj + 1e-300),
         direction = ifelse(avg_log2FC > 0, "Up", "Down"))

# Genes labeled in the figure (genes absent from the DEG table are not drawn)
genes_to_label <- c(
  "ISL1", "ST18", "SOX2", "NR2F1", "MAFB", "HOXB9", "MNX1", "SOX4",
  "ONECUT2", "ROBO2", "CNTN2", "CDH7", "DLL3", "SEMA3C", "SEMA6A",
  "DLGAP2", "PCDH17", "CLSTN2", "SST", "CHRNA3", "PCSK2", "CACNB2",
  "PIEZO2", "ANGPT1", "THSD4", "MDK", "FNDC5", "METRN", "ADAM12",
  "HMGA2", "HMGB2", "CKS2", "HIST1H4C", "ATP5ME", "NDUFA3", "PRDX1",
  "HMGN2", "POLR2L", "DSCAM", "EFNB2", "MDGA2", "LRRTM2", "NRG1",
  "PCDH10", "PCDH15", "CNTN1", "PCDH9", "NCAM2", "NRG3", "MAPT",
  "SV2C", "DNM1", "SYT14", "APBA1", "DNM3", "CPLX2", "SYP", "SYT17",
  "SCN9A", "PLCB1", "GRIK1", "CACNA1C", "SCN3B", "GRIK2", "GABRA2",
  "CAMK2B", "CACNA1E", "NTRK3", "SCN2A", "GRIA4", "PRKCQ", "SMARCA2",
  "SESN3", "SOCS3", "ALDH1A2", "NEUROD1", "PTPRZ1", "PTCH1", "RBP1",
  "HES6", "RAB13", "MEST", "DHRS3", "PRICKLE1"
)
label_df <- volcano_df %>% filter(gene %in% genes_to_label)

p_volcano <- ggplot(volcano_df, aes(x = avg_log2FC, y = neglog10p)) +
  geom_point(color = "black", size = 1.5, alpha = 0.5) +
  geom_point(data = label_df, aes(color = direction), size = 2.5) +
  geom_text_repel(data = label_df, aes(label = gene, color = direction),
                  size = 3, max.overlaps = Inf, box.padding = 0.4,
                  point.padding = 0.3, segment.size = 0.3,
                  segment.alpha = 0.6, show.legend = FALSE) +
  scale_color_manual(values = c(Down = "blue", Up = "red")) +
  labs(x = "log2 fold change (FLUX / Pulsatile)", y = "-log10(adjusted p)") +
  theme_classic()
ggsave(file.path(results_dir, "Fig6b_Volcano_FLUX_vs_Pulsatile.pdf"),
       p_volcano, width = 12, height = 12)

save_session_info("02_Fig6b")
