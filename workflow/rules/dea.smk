
# perform differential expression analysis
rule dea:
    input:
        get_data_path
    output:
        dea_results = os.path.join(result_path,'{analysis}','results.csv'),
        all_features = os.path.join(result_path,'{analysis}','feature_lists',"ALL_features.txt"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","dea_{analysis}.log"),
    params:
        partition=config.get("partition"),
        assay = lambda w: annot_dict["{}".format(w.analysis)]["assay"],
        metadata = lambda w: annot_dict["{}".format(w.analysis)]["metadata"],
        control = lambda w: annot_dict["{}".format(w.analysis)]["control"],
    script:
        "../scripts/dea.R"

# aggregate results per analysis
rule aggregate:
    input:
        dea_results = os.path.join(result_path,'{analysis}','results.csv'),
    output:
        dea_stats = report(os.path.join(result_path,'{analysis}','stats.csv'), 
                                  caption="../report/dea_stats.rst", 
                                  category="{}_{}".format(config["project_name"], module_name),
                                  subcategory="{analysis}",
                                  labels={
                                      "name": "Statistics",
                                      "type": "table",
                                      "misc": "CSV",
                                  }),
#         dea_filtered_stats = report(os.path.join(result_path,'{analysis}','DEA_FILTERED_stats.csv'), 
#                                   caption="../report/dea_stats.rst", 
#                                   category="{}_dea_seurat".format(config["project_name"]), 
#                                   subcategory="{analysis}"),
#         dea_filtered_lfc = os.path.join(result_path,'{analysis}','DEA_FILTERED_LFC.csv'),
        dea_stats_plot = report(os.path.join(result_path,'{analysis}','plots','stats.png'), 
                                  caption="../report/dea_stats.rst", 
                                  category="{}_{}".format(config["project_name"], module_name),
                                  subcategory="{analysis}",
                                  labels={
                                      "name": "Statistics",
                                      "type": "stacked bar plot",
                                      "misc": "PNG",
                                  }),
#         dea_filtered_stats_plot = report(os.path.join(result_path,'{analysis}','plots','DEA_FILTERED_stats.png'), 
#                                   caption="../report/dea_stats.rst", 
#                                   category="{}_dea_seurat".format(config["project_name"]), 
#                                   subcategory="{analysis}"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","aggregate_{analysis}.log"),
    params:
        partition=config.get("partition"),
        assay = lambda w: annot_dict["{}".format(w.analysis)]["assay"],
        metadata = lambda w: annot_dict["{}".format(w.analysis)]["metadata"],
        control = lambda w: annot_dict["{}".format(w.analysis)]["control"],
    script:
        "../scripts/aggregate.R"

