#### load libraries & utility function 
library(pheatmap)
library(patchwork)
library(ggplot2)
library(ggplotify)

# source utility functions
source("workflow/scripts/utils.R")

# inputs
dea_filtered_lfc_path <- snakemake@input[["dea_filtered_lfc"]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/dea_seurat/KOcall_NonTargeting_condition/DEA_FILTERED_LFC.csv"

# outputs
dea_lfc_heatmap_path <- snakemake@output[["dea_lfc_heatmap"]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/dea_seurat/KOcall_NonTargeting_condition/plots/DEA_LFC_heatmap.png"

# parameters
assay <- snakemake@params[["assay"]] #"SCT" #"RNA"
metadata <- snakemake@params[["metadata"]] #"condition"
control <- snakemake@params[["control"]] #"untreated"


# plot specifications
width <- 0.25
height <- 5


### load LFC DEA results
dea_lfc <- read.csv(file=file.path(dea_filtered_lfc_path), row.names = 1)

# set NA values to 0 (NA because below LFC threshold during testing or filtering)
dea_lfc[is.na(dea_lfc)] <- 0

### visualize LFC of DEA results as heatmap
width_panel <- width * ncol(dea_lfc) + 3

annot <- data.frame(group=gsub("group_","",colnames(dea_lfc)))
rownames(annot) <- colnames(dea_lfc)

# make heatmap
lfc_heatmap <- as.ggplot(pheatmap(dea_lfc, 
               show_rownames=F, 
               show_colnames=F,
               treeheight_row = 25,
               treeheight_col = 10,
               annotation_col = annot,
               breaks=seq(-max(abs(dea_lfc)), max(abs(dea_lfc)), length.out=200),
               color=colorRampPalette(c("blue", "white", "red"))(200),
                                     annotation_names_col = F
              ))


# save plot
options(repr.plot.width=width_panel, repr.plot.height=height)
print(lfc_heatmap)


ggsave_new(filename = "DEA_LFC_heatmap", 
           results_path=dirname(dea_lfc_heatmap_path), 
           plot=lfc_heatmap, 
           width=width_panel, 
           height=height)
