
# perform differential expression analysis
rule dea:
    input:
        get_data_path
    output:
        dea_results = os.path.join(result_path,'{analysis}','results.csv'),
        all_features = os.path.join(result_path,'{analysis}','feature_lists',"ALL_features.txt"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: 8*config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","dea_{analysis}.log"),
    params:
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
        dea_stats_plot = report(os.path.join(result_path,'{analysis}','plots','stats.png'), 
                                  caption="../report/dea_stats.rst", 
                                  category="{}_{}".format(config["project_name"], module_name),
                                  subcategory="{analysis}",
                                  labels={
                                      "name": "Statistics",
                                      "type": "Bar plot",
                                      "misc": "PNG",
                                  }),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","aggregate_{analysis}.log"),
    params:
        group = "ALL",
        assay = lambda w: annot_dict["{}".format(w.analysis)]["assay"],
        metadata = lambda w: annot_dict["{}".format(w.analysis)]["metadata"],
        control = lambda w: annot_dict["{}".format(w.analysis)]["control"],
    script:
        "../scripts/aggregate.R"

# generate feature lists per group
# not part of target rule all and only used when the outputs are required by a subsequent module e.g., enrichment_analysis
rule feature_lists:
    input:
        dea_results = os.path.join(result_path,'{analysis}','results.csv'),
        dea_stats = os.path.join(result_path,'{analysis}','stats.csv'), # this ensures rule order (after rule aggregate) to avoid file writing clashes
    output:
        features_up = os.path.join(result_path,'{analysis}','feature_lists',"{group}_up_features.txt"),
        features_down = os.path.join(result_path,'{analysis}','feature_lists',"{group}_down_features.txt"),
        features_scores = os.path.join(result_path,'{analysis}','feature_lists',"{group}_featureScores.csv"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/ggplot.yaml"
    log:
        os.path.join("logs","rules","feature_lists_{analysis}_{group}.log"),
    params:
        group = lambda w: "{}".format(w.group),
        assay = lambda w: annot_dict["{}".format(w.analysis)]["assay"],
        metadata = lambda w: annot_dict["{}".format(w.analysis)]["metadata"],
        control = lambda w: annot_dict["{}".format(w.analysis)]["control"],
    script:
        "../scripts/aggregate.R"