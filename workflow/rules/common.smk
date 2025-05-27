##### utility functions #####

def get_data_path(wildcards):
    return annot.loc[wildcards.analysis,'data']

def get_feature_list_path(wildcards):
    if wildcards.feature_list != "ALL" and wildcards.feature_list != "FILTERED":
        return config["feature_lists"][wildcards.feature_list]
    else:
        return []