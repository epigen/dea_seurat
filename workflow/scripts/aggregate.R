#### load libraries & utility function 
library("ggplot2")

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# inputs
dea_result_path <- snakemake@input[["dea_results"]]

# outputs
dea_stats_path <- snakemake@output[["dea_stats"]]
dea_stats_plot_path <- snakemake@output[["dea_stats_plot"]]

# parameters
assay <- snakemake@params[["assay"]]
metadata <- snakemake@params[["metadata"]]
control <- snakemake@params[["control"]]

adj_pval <- snakemake@config[["filters"]][["adj_pval"]]
lfc <- snakemake@config[["filters"]][["lfc"]]
min_pct <- snakemake@config[["filters"]][["min_pct"]]
score_formula <- snakemake@config[["score_formula"]]

# plot specifications
width <- 0.25
height <- 4

### load DEA results
# dea_results <- read.csv(file=file.path(dea_result_path))
dea_results <- data.frame(fread(file.path(dea_result_path), header=TRUE))
groups <- unique(dea_results$group)

# determine and save feature scores for each gene and group for downstream analysis e.g., preranked GSEA
if (score_formula!=""){
    dea_results$score <- eval(parse(text=score_formula))
    
    for (group in groups){
        tmp_features <- dea_results[(dea_results$group==group),c("feature","score")]
#         write.csv(tmp_features, file.path(dirname(dea_result_path),paste0(group,"_featureScores.csv")), row.names=FALSE)
        fwrite(as.data.frame(tmp_features), file=file.path(dirname(dea_result_path),"feature_lists",paste0(group,"_featureScores.csv")), row.names=FALSE)
    }
}

# annotate differential direction (up or down)
dea_results$direction <- as.character(lapply(dea_results$avg_log2FC, function(x) if(x>0){"up"}else{"down"}))

# ### aggregate & save DEA statistics by stat. sign.
# tmp_dea_results <- dea_results[dea_results$p_val_adj <= adj_pval, ]
# dea_stats <- table(tmp_dea_results$group, tmp_dea_results$direction)
# dea_stats_df <- as.data.frame.matrix(dea_stats)
# dea_stats_df$total <- rowSums(dea_stats_df)

# # write.csv(dea_stats_df, file=file.path(dea_stats_path), row.names=TRUE)
# fwrite(as.data.frame(dea_stats_df), file=file.path(dea_stats_path), row.names=TRUE)

### aggregate & save FILTERED DEA statistics
dea_filtered_results <- dea_results[(dea_results$p_val_adj<=adj_pval) & 
                                  (abs(dea_results$avg_log2FC)>=lfc) & 
                                  (apply(dea_results[,c("pct.1","pct.2")], 1, max) >= min_pct), ]
                                                                                          
dea_filtered_stats <- table(dea_filtered_results$group, dea_filtered_results$direction)
dea_filtered_stats_df <- as.data.frame.matrix(dea_filtered_stats)
dea_filtered_stats_df$total <- rowSums(dea_filtered_stats_df)
                                             
# write.csv(dea_filtered_stats_df, file=file.path(dea_stats_path), row.names=TRUE)
fwrite(as.data.frame(dea_filtered_stats_df), file=file.path(dea_stats_path), row.names=TRUE)

### aggregate & save LFC matrix from filtered DEA results -> NOT USED anymore as it duplicates results unnecessarily
# groups <- paste0("group_",unique(dea_filtered_results$group))
# features <- unique(dea_filtered_results$feature)
                                             
# lfc_df <- data.frame(matrix(nrow=length(features), ncol=length(groups), dimnames=list(features, groups)))

# for (group in unique(dea_filtered_results$group)){
#     tmp_dea_results <- dea_results[dea_results$group==group, ] #old: dea_filtered_results[dea_filtered_results$group==group, ]
#     rownames(tmp_dea_results) <- tmp_dea_results$feature
#     lfc_df[features, paste0("group_",group)] <- tmp_dea_results[features, 'avg_log2FC']
# }

# # write.csv(lfc_df, file=file.path(dea_filtered_lfc_path), row.names=TRUE)
# fwrite(as.data.frame(lfc_df), file=file.path(dea_filtered_lfc_path), row.names=TRUE)

### save differential feature lists from filtered DEA results for downstream analysis (e.g., enrichment analysis)         
for (group in unique(dea_filtered_results$group)){
    for (direction in unique(dea_filtered_results$direction)){
        
        tmp_features <- dea_filtered_results[(dea_filtered_results$group==group) & (dea_filtered_results$direction==direction),"feature"]
        write(tmp_features, file.path(dirname(dea_result_path),"feature_lists",paste0(group,"_",direction,"_features.txt")))
    }
}
                                             
# ### visualize & save ALL DEA statistics -> NOT USED anymore as it's confusing/redundant
# width_panel <- length(groups) * width + 1
                                             
# # convert group column to string for correct plotting
# dea_results$group <- as.character(dea_results$group)
# tmp_dea_results <- dea_results[dea_results$p_val_adj<=adj_pval, ]

# dea_results_p <- ggplot(tmp_dea_results, aes(x=group, fill=direction)) + 
#                                              geom_bar() + 
#                                              xlab(metadata) +
#                                              ylab("number of differential features") +
#                                              scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
#                                              scale_fill_manual(values=list(up="red", down="blue"), drop=FALSE) +
#                                              custom_theme +
#                                              theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size = 6)) # rotates the x-Axis
                                             
                                             
# # save plot
# # options(repr.plot.width=width_panel, repr.plot.height=height)
# # print(dea_results_p)
# ggsave_new(filename = "stats", 
#            results_path=dirname(dea_stats_plot_path), 
#            plot=dea_results_p, 
#            width=width_panel, 
#            height=height)      
                                             
### visualize & save FILTERED DEA statistics
width_panel <- length(groups) * width + 1

# convert group column to string for correct plotting
# dea_filtered_results$group <- as.character(dea_filtered_results$group)
                                             
# dea_filtered_results_p <- ggplot(dea_filtered_results, aes(x=group, fill=direction)) + 
#                                              geom_bar() + 
#                                              xlab(metadata) +
#                                              ylab("number of differential features") +
#                                              scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
#                                              scale_fill_manual(values=list(up="red", down="blue"), drop=FALSE) +
#                                              custom_theme +
#                                              theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size = 6)) # rotates the x-Axis

# fromat stats df for plotting: remove "total" column & transform from wide to long
dea_filtered_stats_df <- dea_filtered_stats_df[,-3]
dea_filtered_stats_df[,"down"] <- -1 * dea_filtered_stats_df[,"down"]
plot_stats_df <- stack(dea_filtered_stats_df)
colnames(plot_stats_df) <- c("n_features","direction")
plot_stats_df$groups <- c(rownames(dea_filtered_stats_df), rownames(dea_filtered_stats_df))

# plot
dea_filtered_results_p <- ggplot(plot_stats_df, aes(x=groups, y=n_features, fill=direction)) + 
                                             geom_bar(stat="identity", position="identity") +
                                             xlab(metadata) +
                                             ylab("number of differential features") +
                                             scale_fill_manual(values=list("down"="blue", "up"="red"), drop=FALSE) +
                                             scale_y_continuous(labels = function(y) sapply(y, function(y) ifelse(y < 0, paste0(sub("-", "", as.character(y))), y))) +
                                             custom_theme +
                                             theme(#legend.position = "none",
                                                   axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6)
                                                  )
                                             
                                             
# save plot
# options(repr.plot.width=width_panel, repr.plot.height=height)
# print(dea_filtered_results_p)
ggsave_new(filename = "stats", 
           results_path=dirname(dea_stats_plot_path), 
           plot=dea_filtered_results_p, 
           width=width_panel, 
           height=height) 
