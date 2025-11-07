# WGCNA Microarray Analysis Pipeline

This repository contains an R-based workflow for performing **Weighted Gene Co-expression Network Analysis (WGCNA)** on microarray expression data.

# Overview
The project includes:
- Normalization and QC steps
- Module detection using WGCNA
- Module–trait correlation analysis
- Generation of `geneInfo.csv` and module–trait heatmap (`moduleTraitHeatmap.png`)

# Files
- `wgcna_pipeline.R` — Main R script containing all analysis steps
- `normalized_expression_matrix.csv` — Expression data input (not uploaded due to size)
- `trait_data.xlsx` — Trait metadata used for correlation (not uploaded)
- `results/` — Folder for heatmaps and output files

# How to Run
1. Clone this repo or download it as ZIP.
2. Place your `normalized_expression_matrix.csv` and `trait_data.xlsx` files in your working directory.
3. Open `wgcna_pipeline.R` in RStudio.
4. Run the script line by line.

# Results
The results of the WGCNA analysis are stored in the `results/` folder:
- `moduleTraitHeatmap.png` — correlation between gene modules and traits
- `geneInfo.csv` — list of genes with module assignments and significance
- `dendro_colors.png` — hierarchical clustering and module colors
- `softThreshold_plot.png` — power selection diagnostics

# Requirements
- R version ≥ 4.2
- Packages: `WGCNA`, `readxl`, `tidyverse`

# Author
Shreemaran V  
📧 shreemaranv13@gmail.com  
🔗 GitHub: shreemaranv-hub(https://github.com/shreemaranv-hub)
