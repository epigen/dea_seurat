##### utility functions #####

def get_data_path(wildcards):
    return annot.loc[wildcards.analysis,'data']

def get_feature_list_path(wildcards):
    return config["feature_lists"][wildcards.feature_list]