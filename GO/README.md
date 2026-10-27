# GO/

Gene Ontology (GO) biological process terms underlying Fig. 6c-f (J. Seong). These tables are identical to Supplementary Data 1-4 of the article.

| File | Content | Figure |
|---|---|---|
| Supplementary Data 1 | GO terms for genes upregulated in FLUX MNs, grouped into 15 functional categories (fold enrichment > 10%) | 6c |
| Supplementary Data 2 | GO terms for genes upregulated in Pulsatile MNs, grouped into 15 functional categories (fold enrichment > 10%) | 6c |
| Supplementary Data 3 | FLUX GO terms relevant to the nervous system and MN function, assigned to four developmental stages | 6d-f |
| Supplementary Data 4 | Pulsatile GO terms relevant to the nervous system and MN function, assigned to four developmental stages | 6d-f |

Input gene lists are the DEGs produced by `R/02_Fig6b_DEG_FLUX_vs_Pulsatile.R`. GO enrichment was run with clusterProfiler `enrichGO` (biological process); grouping into functional categories and developmental stages was done manually as described in the Methods.
