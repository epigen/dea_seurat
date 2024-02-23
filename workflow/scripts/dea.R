#### load libraries & utility function 
library("Seurat")
library("data.table")
library("future")

# change the current plan to access all available cores for parallelization
plan("multicore") # "multisession" does not work
plan()
# set global size for parallelization -> TODO: dynamic depending on data size, RAM and available threads. how?
options(future.globals.maxSize = 1000 * 1024^2)

# inputs
object_path <- snakemake@input[[1]]

# outputs
dea_result_path <- snakemake@output[["dea_results"]]
all_features_path <- snakemake@output[["all_features"]]

# parameters
assay <- snakemake@params[["assay"]]
metadata <- snakemake@params[["metadata"]]
control <- snakemake@params[["control"]]

logfc_threshold <- snakemake@config[["logfc_threshold"]]
test_use <- snakemake@config[["test_use"]]
min_pct <- snakemake@config[["min_pct"]]
return_thresh <- snakemake@config[["return_thresh"]]

### load data
data <- readRDS(file = file.path(object_path))
DefaultAssay(object = data) <- assay
Idents(object = data) <- metadata

# Prepare object to run differential expression on SCT assay with multiple models (just in case)
if(assay=="SCT"){
    data <- PrepSCTFindMarkers(data, assay = assay, verbose = TRUE)
}

### save list of all expressed features for downstream analysis (e.g., as background in enrichment analyses)
all_features <- rownames(GetAssayData(object = data, assay = assay, slot = "data"))
write(all_features, file.path(all_features_path))

### perform DEA
if (control=="ALL"){
    # one vs all DEA to identify group specific markers compared to rest
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
    # group vs control DEA to identify changes compared to control
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
fwrite(as.data.frame(dea_results), file=file.path(dea_result_path), row.names=FALSE)