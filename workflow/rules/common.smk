##### utility functions #####

def get_data_path(wildcards):
    return annot.loc[wildcards.analysis,'data']
