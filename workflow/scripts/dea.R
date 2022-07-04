#### load libraries & utility function 
library(Seurat)

# inputs
object_path <- snakemake@input[[1]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/scrnaseq_processing_seurat/KOcall_NonTargeting/NORMALIZED_object.rds"

# outputs
dea_result_path <- snakemake@output[["dea_results"]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/dea_seurat/KOcall_NonTargeting_condition/DEA_results.csv" 

# parameters
assay <- snakemake@params[["assay"]] #"SCT" #"RNA"
metadata <- snakemake@params[["metadata"]] #"condition"
control <- snakemake@params[["control"]] #"untreated"

logfc_threshold <- snakemake@params[["logfc_threshold"]] #0.1
test_use <- snakemake@params[["test_use"]] #"wilcox"
min_pct <- snakemake@params[["min_pct"]] #0.1

return_thresh <- snakemake@params[["return_thresh"]] #1


result_dir <- dirname(dea_result_path)
# make directories if not exist
if (!dir.exists(result_dir)){
        dir.create(result_dir, recursive = TRUE)
    }

### load data
data <- readRDS(file = file.path(object_path))
DefaultAssay(object = data) <- assay
Idents(object = data) <- metadata


### save list of all expressed features for downstream analysis (eg as background in enrichment analyses)
features_path <- file.path(result_dir,"feature_lists")
if (!dir.exists(features_path)){
    dir.create(features_path, recursive = TRUE)
}
all_features <- rownames(GetAssayData(object = data, assay = assay, slot = "data"))
write(all_features, file.path(features_path, "ALL_features.txt"))


### perform DEA

# one vs all DEA to identify group specific markers compared to rest
if (control=="ALL"){
    dea_results <- FindAllMarkers(object = data,
                                  assay = assay,
                                  features = NULL,
                                  logfc.threshold = logfc_threshold,
                                  test.use = test_use,
                                  slot = "data",
                                  min.pct = min_pct,
                                  min.diff.pct = -Inf,
                                  node = NULL,
                                  verbose = TRUE,
                                  only.pos = FALSE,
                                  max.cells.per.ident = Inf,
                                  random.seed = 42,
                                  latent.vars = NULL,
                                  min.cells.feature = 3,
                                  min.cells.group = 3,
                                  pseudocount.use = 1,
                                  mean.fxn = NULL,
                                  fc.name = NULL,
                                  base = 2,
                                  return.thresh = return_thresh,
                                  densify = FALSE
                                 )

    
    # rename columns and rows, if results are not empty
    if (dim(dea_results)[1]!=0){
        colnames(dea_results)[colnames(dea_results) == "gene"] <- "feature"
        colnames(dea_results)[colnames(dea_results) == "cluster"] <- "group"
        rownames(dea_results) <- NULL
    }
    
    
}else{
    # grouo vs control DEA to identify changes compared to control
    dea_results <- data.frame()

    for (group in unlist(unique(data[[metadata]]))){
        if (group==control){
            next
        }

        tmp_markers <- FindMarkers(object = data,
                                   ident.1 = group,
                                   ident.2 = control,
                                   group.by = NULL,
                                   subset.ident = NULL,
                                   assay = assay,
                                   slot = "data",
                                   reduction = NULL,
                                   features = NULL,
                                   logfc.threshold = logfc_threshold,
                                   test.use = test_use,
                                   min.pct = min_pct,
                                   min.diff.pct = -Inf,
                                   verbose = TRUE,
                                   only.pos = FALSE,
                                   max.cells.per.ident = Inf,
                                   random.seed = 42,
                                   latent.vars = NULL,
                                   min.cells.feature = 3,
                                   min.cells.group = 3,
                                   pseudocount.use = 1,
                                   mean.fxn = NULL,
                                   fc.name = NULL,
                                   base = 2,
                                   densify = FALSE
                                  )

        # if results are empty move on
        if(dim(tmp_markers)[1]==0){
            next
        }
        
        tmp_markers$group <- group
        tmp_markers$feature <- rownames(tmp_markers)
        rownames(tmp_markers) <- NULL

        if(dim(dea_results)[1]==0){
            dea_results <- tmp_markers
        }else{
            dea_results <- rbind(dea_results, tmp_markers)
        }
    }
}

### save results
write.csv(dea_results, file=file.path(dea_result_path), row.names=FALSE)