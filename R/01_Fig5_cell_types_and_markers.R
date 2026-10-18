# =============================================================================
# 01_Fig5_cell_types_and_markers.R
# Figure 5: cell types and motor neuron (MN) markers in day-25 hiPSC-derived
# cultures, Pulsatile vs FLUX.
#   Fig. 5b    UMAP colored by cell type, split by group
#   Fig. 5c    Top 10 marker genes per cell type (heatmap input tables)
#   Fig. 5d    Feature plots of MN markers
#   Fig. 5e    Cell-type composition (% of all cells) per group
#   Fig. 5f    MNs as % of the neuronal population per group
#   Fig. 5g-i  FOXP1, ISL1 and MNX1 (HB9) expression in MNs
# Expected values are noted next to each output (from the manuscript).
# =============================================================================

source("R/00_config.R")

obj_all <- readRDS(path_hipsc_obj)
obj <- subset(obj_all, subset = orig.ident %in% groups_compared)
obj$orig.ident <- factor(as.character(obj$orig.ident), levels = groups_compared)

table(obj$orig.ident)   # expected: Pulsatile 6380, FLUX 8208 cells

# ---- Fig. 5b: UMAP by cell type, split by group -----------------------------
p_umap <- DimPlot(obj, reduction = umap_reduction, group.by = celltype_col,
                  split.by = "orig.ident", label = FALSE) +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(file.path(results_dir, "Fig5b_UMAP_celltype_split_by_group.pdf"),
       p_umap, width = 16, height = 8)

# ---- Fig. 5c: Top 10 marker genes per cell type -----------------------------
# PrepSCTFindMarkers is applied to a copy so that the expression values used
# for Fig. 5g-i come from the unmodified object, as in the original analysis.
obj_markers <- obj
Idents(obj_markers) <- celltype_col
DefaultAssay(obj_markers) <- "SCT"
obj_markers <- PrepSCTFindMarkers(obj_markers, assay = "SCT")

markers_all <- FindAllMarkers(obj_markers, assay = "SCT", slot = "data",
                              only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.25, test.use = "wilcox")
markers_top10 <- markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()
write_csv(markers_top10,
          file.path(results_dir, "Fig5c_top10_markers_by_celltype_SCT.csv"))
rm(obj_markers)

# Average SCT expression per cell type; the Fig. 5c heatmap shows these
# values (scaled per gene) for the top 10 markers. As in the original
# analysis, averages are computed over all three libraries.
obj_all$celltype <- as.character(obj_all[[celltype_col, drop = TRUE]])
obj_ct <- subset(obj_all, subset = !is.na(celltype) & celltype != "")
Idents(obj_ct) <- "celltype"
avg_by_celltype <- AverageExpression(obj_ct, assays = "SCT", slot = "data",
                                     verbose = FALSE)$SCT %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "celltype", values_to = "avg_exp")
write_csv(avg_by_celltype,
          file.path(results_dir, "Fig5c_avg_expression_SCT_by_celltype.csv"))
rm(obj_ct)

# ---- Fig. 5d: Feature plots of canonical MN markers -------------------------
# MNX1 is labeled HB9 in the figure.
mn_features <- c("ISL1", "MNX1", "NEFM", "STMN2", "SLC18A3", "SLC5A7")
feature_plots <- lapply(mn_features, function(gene) {
  FeaturePlot(obj, features = gene, reduction = umap_reduction, pt.size = 0.1) &
    theme(axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank())
})
ggsave(file.path(results_dir, "Fig5d_FeaturePlots_MN_markers.pdf"),
       wrap_plots(feature_plots, ncol = 1), width = 5, height = 24)

# ---- Fig. 5e: Cell-type composition (% of all cells) ------------------------
celltype_percent <- obj@meta.data %>%
  count(orig.ident, .data[[celltype_col]], name = "celltype_count") %>%
  group_by(orig.ident) %>%
  mutate(total_cells = sum(celltype_count),
         celltype_percent = 100 * celltype_count / total_cells) %>%
  ungroup()
write_csv(celltype_percent,
          file.path(results_dir, "Fig5e_celltype_percent_by_group.csv"))
# expected (FLUX vs Pulsatile): Pan-N 14.32 vs 17.96; NPC 56.02 vs 54.37;
# MNP 4.34 vs 3.95; GABA/GlyN 0.67 vs 1.11; SN 0.06 vs 0.13

# ---- Fig. 5f: MNs as % of neuronal cells ------------------------------------
neuron_percent <- celltype_percent %>%
  filter(.data[[celltype_col]] %in% neuron_celltypes) %>%
  group_by(orig.ident) %>%
  mutate(neuron_total = sum(celltype_count),
         percent_of_neurons = 100 * celltype_count / neuron_total) %>%
  ungroup()
write_csv(neuron_percent,
          file.path(results_dir, "Fig5f_percent_of_neuronal_population.csv"))
# expected MN share: FLUX 42.42, Pulsatile 34.73

# MNX1 (HB9)-expressing cells as % of all cells (reported in the text, ~5%)
# Computed on the default assay of the object, as originally run
expr_default <- GetAssayData(obj, assay = DefaultAssay(obj), layer = "data")
hb9_percent <- obj@meta.data %>%
  mutate(hb9_pos = expr_default["MNX1", rownames(obj@meta.data)] > 0) %>%
  group_by(orig.ident) %>%
  summarise(percent_MNX1_expressing = 100 * mean(hb9_pos), .groups = "drop")
write_csv(hb9_percent,
          file.path(results_dir, "Text_percent_MNX1_expressing_cells.csv"))
# expected: FLUX 6.06, Pulsatile 4.50

# ---- Fig. 5g-i: FOXP1, ISL1 and MNX1 in MNs ---------------------------------
marker_genes <- c("FOXP1", "ISL1", "MNX1")

# Mean SCT expression in somatic MNs (values quoted in the text)
obj_smn <- subset(obj, subset = new_consolidated_annotation == somatic_mn)
Idents(obj_smn) <- "orig.ident"
mn_means <- AverageExpression(obj_smn, assays = "SCT", slot = "data",
                              features = marker_genes, verbose = FALSE)$SCT %>%
  as.data.frame() %>%
  rownames_to_column("gene")
write_csv(mn_means, file.path(results_dir, "Fig5g-i_mean_SCT_expression_MN.csv"))
# expected (FLUX vs Pulsatile): FOXP1 3.74 vs 3.37; ISL1 0.93 vs 0.44;
# MNX1 0.64 vs 0.46

# Violin plots (cells from the MN and v-MN clusters, as in the figure file)
obj_mn <- subset(obj, subset = new_consolidated_annotation %in% mn_celltypes)
p_violin <- VlnPlot(obj_mn, features = marker_genes, group.by = "orig.ident",
                    assay = "SCT", layer = "data", pt.size = 0.1, ncol = 3)
ggsave(file.path(results_dir, "Fig5g-i_Violin_FOXP1_ISL1_MNX1.pdf"),
       p_violin, width = 9, height = 6)

save_session_info("01_Fig5")
