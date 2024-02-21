# Configuration

You need one configuration file and one annotation file to run the complete workflow. You can use the provided example as starting point. If in doubt read the comments in the config and/or try the default values.

- project configuration (config/config.yaml): different for every project/dataset and configures the analyses to be performed.
- sample annotation (sample_annotation): CSV file consisting of five columns
    -  name: name of the dataset/analysis (tip: keep it short, but descriptive and distinctive).
    -  data: path to the input Seurat object as .rds.
    -  assay: the Seurat assay to be used (e.g., SCT or RNA).
    -  metadata: column name of the metadata that should be used to group cells for comparison (e.g., condition or cell_type).
    -  control: name of the class/level that should be used as control in the comparison (e.g., untreated) or "ALL" to compare every class against the rest (e.g., useful to find cluster markers; one vs all)
