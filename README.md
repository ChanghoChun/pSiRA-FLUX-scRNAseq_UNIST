# FLUX scRNA-seq analysis code

Code for the single-cell RNA-seq analyses in:

> Seong J.\*, Chun C.\*, Kim M., Bae Y., Park J., Kim J., Kang J.H., Park T.-E., Bothwell M., Smith A.S.T., Mack D.L., Joo J. FLUX: Facilitated long-term uninterrupted exposure of retinoic acid enhances transcriptomic and electrophysiological maturation of iPSC-derived motor neurons. *Nature Communications* (accepted in principle). \*Equal contribution.

To cite this code, use the Zenodo DOI: **[Zenodo DOI]** (see `CITATION.cff`).

## Data

| Data | Source |
|---|---|
| Raw FASTQ files and Cell Ranger filtered matrices (3 libraries, day 25) | NCBI GEO, accession **[GSE number]** |
| Human embryonic spinal cord scRNA-seq (Rayon et al., 2021) | GEO GSE171892 |
| Second human spinal cord reference for scANVI (step 2) | [to be added with `python/`] |

Library names used in the code (`orig.ident`):

| Code | Condition | GEO sample | FASTQ prefix |
|---|---|---|---|
| `Pulsatile` | Pulsatile RA | RA_only | Ctrl-MN |
| `FLUX` | Sustained RA release from RA-loaded pSiMPs | RA_at_pSiMP | RA-MN |
| `pSiMP` | Pulsatile RA + blank pSiMPs | RA_plus_pSiMP | pSi-MN |

## Workflow

1. **Alignment**: Cell Ranger v7.0.1, reference GRCh38-2020-A. FASTQ to filtered gene-barcode matrices (deposited in GEO).
2. **Quality control, integration and annotation** (`python/`, to be added): doublet removal (scDblFinder), MAD-based filtering (scuttle), ambient RNA correction (SoupX), SCTransform normalization, scANVI reference model and label transfer (scvi-tools v1.1.5). Output: the two Seurat objects listed in `R/00_config.R`, placed in `data/`.
3. **Figure analyses** (`R/`):

| Script | Figure panels | Main outputs |
|---|---|---|
| `00_config.R` | | Paths, packages, shared definitions |
| `01_Fig5_cell_types_and_markers.R` | 5b, 5c, 5d, 5e, 5f, 5g-i | UMAP, marker tables, feature plots, composition tables, MN marker means and violin plots |
| `02_Fig6b_DEG_FLUX_vs_Pulsatile.R` | 6b | DEG table (FLUX vs Pulsatile MNs), volcano plot |
| `03_Fig8_Spearman_vs_human_embryonic_MN.R` | 8b-g | Spearman correlation table (hiPSC vs embryonic somatic MNs) |
| `GO/` (to be added) | 6c-f, 7a-d | GO enrichment and category-level analyses |

Each script lists the values reported in the manuscript next to the corresponding output, so results can be checked directly.

Run from the repository root:

```bash
Rscript R/01_Fig5_cell_types_and_markers.R
Rscript R/02_Fig6b_DEG_FLUX_vs_Pulsatile.R
Rscript R/03_Fig8_Spearman_vs_human_embryonic_MN.R
```

Results are written to `results/`. Bar charts and heatmaps in the final figures were drawn in GraphPad Prism from the exported tables. Schematic panels (5a, 6a, 6d, 7e, 8a) contain no data.

## Software

R 4.4.2 with Seurat 5.3.0, SeuratObject 5.1.99.9000, dplyr 1.1.4, tidyr 1.3.1, tibble 3.3.0, readr 2.1.5, ggplot2 3.5.2, ggrepel 0.9.6 and patchwork 1.3.1. All figure values listed in the scripts were reproduced with these versions (macOS). Python 3.11, scvi-tools 1.1.5, scanpy 1.10.4 and anndata 0.11.3 were used for step 2. Each script saves the full session information to `results/sessionInfo_*.txt`.

## License

MIT License. See `LICENSE`.
