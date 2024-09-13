
# visualize DEA using volcano plots
rule volcanos:
    input:
        results = os.path.join(result_path,'{analysis}','results.csv'),
    output:
        volcanos = report(directory(os.path.join(result_path,'{analysis}','plots','volcano','{feature_list}')),
                              patterns=["{group}.png"],
                              caption="../report/volcano.rst",
                              category="{}_{}".format(config["project_name"], module_name),
                              subcategory="{analysis}",
                              labels={
                                      "name": "Volcano",
                                      "type": "{group}",
                                      "misc": "{feature_list}",
                                  }),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/volcanos.yaml"
    log:
        os.path.join("logs","rules","volcanos_{analysis}_{feature_list}.log"),
    params:
        assay = lambda w: annot_dict["{}".format(w.analysis)]["assay"],
        metadata = lambda w: annot_dict["{}".format(w.analysis)]["metadata"],
        control = lambda w: annot_dict["{}".format(w.analysis)]["control"],
    script:
        "../scripts/volcanos.R"
        
# visualize LFC of DEA results
rule heatmap:
    input:
        results = os.path.join(result_path,'{analysis}','results.csv'),
    output:
        lfc_heatmap = report(os.path.join(result_path,'{analysis}','plots','heatmap','{feature_list}.png'),
                              caption="../report/heatmap.rst",
                              category="{}_{}".format(config["project_name"], module_name),
                              subcategory="{analysis}",
                              labels={
                                      "name": "Heatmap",
                                      "type": "Effect sizes",
                                      "misc": "{feature_list}",
                                  }),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/heatmap.yaml"
    log:
        os.path.join("logs","rules","lfc_heatmap_{analysis}_{feature_list}.log"),
    params:
        assay = lambda w: annot_dict["{}".format(w.analysis)]["assay"],
        metadata = lambda w: annot_dict["{}".format(w.analysis)]["metadata"],
        control = lambda w: annot_dict["{}".format(w.analysis)]["control"],
    script:
        "../scripts/heatmap.R"