#### load libraries & utility function 
library("EnhancedVolcano", quietly=TRUE)
library("patchwork", quietly=TRUE)
library("ggplot2")

options(ragg.max_dim = 100000) # required for large volcano panels

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# inputs
dea_result_path <- snakemake@input[["results"]]

# outputs
volcano_plot_path <- snakemake@output[["volcanos"]]

# parameters
assay <- snakemake@params[["assay"]]
metadata <- snakemake@params[["metadata"]]
control <- snakemake@params[["control"]]

pCutoff <- snakemake@config[["volcano"]][["pCutoff"]]
FCcutoff <- snakemake@config[["volcano"]][["FCcutoff"]]
feature_list_name <- snakemake@wildcards[["feature_list"]]


# plot specifications
# n_col <- 10
width <- 4
height <- 5
# width_panel <- n_col * width


### load DEA results
# dea_results <- read.csv(file=file.path(dea_result_path))
dea_results <- data.frame(fread(file.path(dea_result_path), header=TRUE))

# load feature list if not ALL
if (feature_list_name!="ALL"){
    feature_list_path <- snakemake@config[["feature_lists"]][[feature_list_name]]
    feature_list <- scan(file.path(feature_list_path), character())
}

# convert group column to string for correct plotting
dea_results$group <- as.character(dea_results$group)

# height_panel <- height * ceiling(length(unique(dea_results$group))/n_col)

### Visualize DEA results using Volcano plots
# volcano_plots <- list()
for (pval_type in c("p_val_adj", "p_val")){
    for (group in unique(dea_results$group)){
        
        toptable <- dea_results[dea_results$group==group,]
        lab <- toptable$feature
        x <- "avg_log2FC"
        
        # set (adjusted) p-values of 0 with min(adj.P.Val>0) * 10ˆ-1 for plotting
        toptable[toptable$p_val_adj==0, "p_val_adj"] <- min(toptable$p_val_adj[toptable$p_val_adj!=0])*10^-1
        toptable[toptable$p_val==0, "p_val"] <- min(toptable$p_val[toptable$p_val!=0])*10^-1

        volcano_plot <- EnhancedVolcano(toptable = toptable,
                        lab = lab,
                        x = x,
                        y = pval_type,
                        selectLab = NULL,
                        xlim = c(min(toptable[[x]], na.rm = TRUE) - 1, max(toptable[[x]], na.rm = TRUE) + 1),
                        ylim = c(0, max(-log10(toptable[[y]]), na.rm = TRUE) + 5),
                        xlab = bquote("average" ~log[2] ~ "fold change"),
                        ylab = if(pval_type=='p_val_adj') bquote(~-log[10] ~ "adjusted p-value") else bquote(~-log[10] ~ "raw p-value"),
                        axisLabSize = 12,
                        title = paste0(metadata,': ', group,' vs ', control),
                        subtitle = paste0('assay: ',assay), # default: bquote(italic(EnhancedVolcano))
                        caption = paste0("variables:",nrow(toptable),"; |avg.log2FC|>",FCcutoff,"; adj.p-val<",pCutoff),
                        titleLabSize = 14,
                        subtitleLabSize = 14,
                        captionLabSize = 6,
                        pCutoff = pCutoff, #default: 0.05
                        pCutoffCol = "p_val_adj", # cutoff always using adjusted p-values
                        FCcutoff = FCcutoff, # default:1
                        cutoffLineType = "longdash",
                        cutoffLineCol = "black",
                        cutoffLineWidth = 0.4,
                        pointSize = 1, # default: 2
                        labSize = 3, #default: 5
                        labCol = "black",
                        labFace = "plain",
                        boxedLabels = TRUE, #default: FALSE
                        parseLabels = FALSE,
                        shape = 19,
                        shapeCustom = NULL,
                        col = c("grey30", "forestgreen", "royalblue", "red2"),
                        colCustom = NULL,
                        colAlpha = 1/2,
                        colGradient = NULL,
                        colGradientBreaks = c(pCutoff, 1),
                        colGradientLabels = c("0", "1.0"),
                        colGradientLimits = c(0, 1),
                        legendLabels = c("NS", expression(avg. ~ log[2] ~ FC), "adj. p-value", 'both'),
                        legendPosition = "right", #default: "top"
                        legendLabSize = 14,
                        legendIconSize = 4,
                        legendDropLevels = TRUE,
                        encircle = NULL,
                        encircleCol = "black",
                        encircleFill = "pink",
                        encircleAlpha = 3/4,
                        encircleSize = 2.5,
                        shade = NULL,
                        shadeFill = "grey",
                        shadeAlpha = 1/2,
                        shadeSize = 0.01,
                        shadeBins = 2,
                        drawConnectors = TRUE, #default: FALSE
                        widthConnectors = 0.1, # default: 0.5
                        typeConnectors = "closed",
                        endsConnectors = "first",
                        lengthConnectors = unit(0.01, "npc"),
                        colConnectors = "grey10",
                        max.overlaps = 0,
                        maxoverlapsConnectors = 10, # default: NULL
                        min.segment.length = 0,
                        directionConnectors = "both",
                        arrowheads = FALSE, # default: TRUE
                        hline = NULL,
                        hlineType = "longdash",
                        hlineCol = "black",
                        hlineWidth = 0.4,
                        vline = NULL,
                        vlineType = "longdash",
                        vlineCol = "black",
                        vlineWidth = 0.4,
                        gridlines.major = TRUE,
                        gridlines.minor = TRUE,
                        border = "partial",
                        borderWidth = 0.8,
                        borderColour = "black",
                        raster = FALSE
                       ) + custom_theme + theme(legend.title = element_blank())

        set.seed(42)
        ggsave_new(filename = if(pval_type=='p_val_adj') paste0(group,"__adjp") else paste0(group,"__rawp"),
               results_path=volcano_plot_path, 
               plot=volcano_plot, 
               width=width, 
               height=height)

}

# volcano_plots_panel <- wrap_plots(volcano_plots, ncol = n_col, guides = "collect")

# save plot
# options(repr.plot.width=width_panel, repr.plot.height=height_panel)
# print(volcano_plots_panel)

# set.seed(42)
# ggsave_new(filename = "DEA_volcanos", 
#            results_path=dirname(volcano_plot_path), 
#            plot=volcano_plots_panel, 
#            width=width_panel, 
#            height=height_panel)
