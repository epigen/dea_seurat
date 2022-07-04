# scRNA-seq Differential Expression Analysis & Visualization Snakemake Workflow powered by Seurat
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow for performing differential expression analyses (DEA) of (multimodal) sc/snRNA-seq data powered by the R package [Seurat's](https://satijalab.org/seurat/index.html) functions [FindMarkers](https://satijalab.org/seurat/reference/findmarkers) and [FindAllMarkers](https://satijalab.org/seurat/reference/findallmarkers).

**If you use this workflow in a publication, don't forget to give credits to the authors by citing the URL of this (original) repository (and its DOI, see Zenodo badge above -> coming soon).**

![Workflow Rulegraph](./workflow/dags/rulegraph.svg)

Table of contents
----------------
  * [Authors](#authors)
  * [Software](#software)
  * [Methods](#methods)
  * [Features](#features)
  * [Usage](#usage)
  * [Configuration](#configuration)
  * [Examples](#examples)
  * [Links](#links)

# Authors
- [Stephan Reichl](https://github.com/sreichl)


# Software
This project wouldn't be possible without the following software and their dependencies:

| Software       | Reference (DOI)                                   |
| :------------: | :-----------------------------------------------: |
| EnhancedVolcano| https://doi.org/10.18129/B9.bioc.EnhancedVolcano  |
| ggplot2        | https://ggplot2.tidyverse.org/                    |
| patchwork      | https://CRAN.R-project.org/package=patchwork      |
| pheatmap       | https://cran.r-project.org/package=pheatmap       |
| Seurat         | https://doi.org/10.1016/j.cell.2021.04.048        |
| Snakemake      | https://doi.org/10.12688/f1000research.29032.2    |

# Methods
This is a template for the Methods section of a scientific publication and is intended to serve as a starting point. Only retain paragraphs relevant to your analysis. References [ref] to the respective publications are curated in the software table above. Versions (ver) have to be read out from the respective conda environment specifications (.yaml file) or post execution. Parameters that have to be adapted depending on the data or workflow configurations are denoted in squared brackets e.g. [X].

The outlined analyses were performed using the R package Seurat (ver) [ref] unless stated otherwise.

**Differential Expression Analysis (DEA).** DEA was performed on the assay [X] and data slot [X] with Seurat's FindMarkers/FindAllMarkers function using the statistical test [X] with the parameters log2(fold change) threshold of [X] and minimal percentage of expression [X]. The results were filtered for relevant features by adjusted p-value of [X], absolute log2(fold change) of [X] and minimum percentage of expression [X].

**Visualization.** All and filtered result statistics, i.e., number of statistically significant results split by positive (up) and negative (down) effect-sizes, were separately visualized with stacked bar plots using ggplot (ver) [ref]. 
To visually summarize results of the same analysis the filtered log2(fold change) values of features that were found to be at least in one comparison statistically significantly differentially expressed were visualized in a hierarchically clustered heatmap using pheatmap (ver) [ref]. 
Volcano plots were generated for each analysis using EnhancedVolcano (ver) [ref] with adjusted p-value threshold of [X] and log2(fold change) threshold of [X] as visual cut-offs for the y- and x-axis, respectively.

**The analysis and visualizations described here were performed using a publicly available Snakemake [ver] (ref) workflow [ref - cite this workflow here].**

# Features
The workflow perfroms the following steps.
- Differential Expression Analysis (DEA)
  - using Seurat's [FindMarkers](https://satijalab.org/seurat/reference/findmarkers) or [FindAllMarkers](https://satijalab.org/seurat/reference/findallmarkers) depending on the configuration (CSV)
  - feature list per comparison group and direction (up/down) for downstream analysis (eg enrichment analysis) (TXT)
- DEA result statistics: number of statistically significant results split by positive (up) and negative (down) change (CSV)
- DEA result filtering by 
  - statistical significance (adjusted p-value)
  - effect-size (log 2 fold change)
  - expression (minimum percentage of expression) in one of the comparison groups
- Log Fold Change (LFC) matrix of filtered features by comparison groups (CSV)
- Visualizations
  - all and filtered DEA result statistics: number of features and direction (stacked Bar plots)
  - Volanco plot per comparison with configured cutoffs for statistical significance and effect-size
  - Clustered Heatmaps of the LFC matrix

# Usage
Here are some tips for the usage of this workflow:
- perform your first run with loose filtering options/cut-offs and set all the same to see if further filtering is even necessar or useful

# Configuration
Detailed specifications can be found here [./config/README.md](./config/README.md)

# Examples
--- COMING SOON ---

# Links
- [GitHub Repository](https://github.com/epigen/dea_seurat/)
- [GitHub Page](https://epigen.github.io/dea_seurat/)
- [Zenodo Repository (coming soon)]()
- [Snakemake Workflow Catalog Entry](https://snakemake.github.io/snakemake-workflow-catalog?usage=epigen/dea_seurat)
