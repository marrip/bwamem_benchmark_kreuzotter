rule bm_plot:
    input:
        get_bm_files(),
    output:
        report(
            "analysis_output/wes/plot.pdf",
            caption="../report/bm_plot.rst",
            category="Time",
        ),
    log:
        "analysis_output/wes/plot.log",
    container:
        config["tools"]["r"]
    script:
        "../scripts/bm_plot.R"
