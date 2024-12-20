#### load libraries & utility function 
library("pheatmap")
library("ggplotify")
library("reshape2")

# source utility functions
# source("workflow/scripts/utils.R")
# snakemake@source("./utils.R") # does not work when loaded as module (https://github.com/snakemake/snakemake/issues/2205)
source(snakemake@params[["utils_path"]])

# inputs
dea_result_path <- snakemake@input[["results"]]

# outputs
lfc_heatmap_path <- snakemake@output[["lfc_heatmap"]]

# parameters
assay <- snakemake@params[["assay"]]
metadata <- snakemake@params[["metadata"]]
control <- snakemake@params[["control"]]

adj_pval <- snakemake@config[["filters"]][["adj_pval"]]
lfc <- snakemake@config[["filters"]][["lfc"]]
min_pct <- snakemake@config[["filters"]][["min_pct"]]

feature_list_name <- snakemake@wildcards[["feature_list"]]

# plot specifications
width <- 0.15
height <- 0.15

# load dea results
dea_results <- data.frame(fread(file.path(dea_result_path), header=TRUE))

# generate or load feature list
if (feature_list_name=="FILTERED"){
    feature_list <- unique(dea_results[(dea_results$p_val_adj<=adj_pval) & 
                                  (abs(dea_results$avg_log2FC)>=lfc) & 
                                  (apply(dea_results[,c("pct.1","pct.2")], 1, max) >= min_pct), 'feature'])
} else {
    feature_list_path <- snakemake@config[["feature_lists"]][[feature_list_name]]
    feature_list <- scan(file.path(feature_list_path), character())
    # ensure features are in the results
    feature_list <- unique(intersect(feature_list, dea_results$feature))
}

# if no features in the results, end early with empty plot
if(length(feature_list)==0){
    ggsave_new(filename = feature_list_name, 
           results_path=dirname(lfc_heatmap_path), 
           plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Features not found in DEA results.") + theme_void(), 
           width=4, 
           height=1)
    
    quit(save = "no", status = 0)
}

# make LFC dataframe
lfc_df <- reshape2::dcast(dea_results, feature ~ group, value.var = 'avg_log2FC')
rownames(lfc_df) <- lfc_df$feature
lfc_df$feature <- NULL

# make adjusted p-value dataframe
adjp_df <- reshape2::dcast(dea_results, feature ~ group, value.var = 'p_val_adj')
rownames(adjp_df) <- adjp_df$feature
adjp_df$feature <- NULL

# subset dataframes according to feature list
lfc_df <- lfc_df[feature_list, ]
adjp_df <- adjp_df[feature_list, ]

# set NA values to 0 (NA because below LFC threshold during testing or filtering)
lfc_df[is.na(lfc_df)] <- 0

# indicate significance
adjp_df[adjp_df<=adj_pval] <- "*"
adjp_df[adjp_df>adj_pval] <- ""
adjp_df[is.na(adjp_df)] <- ""

### visualize LFC of DEA results as heatmap
height_panel <-  if (nrow(lfc_df)<100) (height * nrow(lfc_df) + 2) else 5
width_panel <- width * ncol(lfc_df) + 2

# make heatmap
lfc_heatmap <- as.ggplot(pheatmap(lfc_df,
                                  display_numbers = if(nrow(lfc_df)<100) adjp_df else FALSE,
                                  main=paste0("Avg.log2FC of ", feature_list_name," features\n", metadata,' vs ', control),
                                  cluster_cols = ifelse(ncol(lfc_df)>1, TRUE, FALSE),
                                  cluster_rows = ifelse(nrow(lfc_df)>1, TRUE, FALSE),
                                  show_rownames = ifelse(nrow(lfc_df)<100, TRUE, FALSE),
                                  labels_row = rownames(lfc_df),
                                  show_colnames = TRUE,
                                  fontsize = 5,
                                  fontsize_number = 10,
                                  angle_col = 45,
                                  treeheight_row = 10,
                                  treeheight_col = 10,
                                  cellwidth = 10,
                                  cellheight = ifelse(nrow(lfc_df)<100, 10, NA),
                                  breaks=seq(-max(abs(lfc_df)), max(abs(lfc_df)), length.out=200),
                                  color=colorRampPalette(c("blue", "white", "red"))(200),
                                  annotation_names_col = F,
                                  silent = TRUE
                                 )
                        )
# save heatmap
ggsave_new(filename = feature_list_name, 
           results_path=dirname(lfc_heatmap_path), 
           plot=lfc_heatmap, 
           width=width_panel, 
           height=height_panel)
