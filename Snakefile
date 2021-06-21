include: "rules/common.smk"
include: "rules/bm_plot.smk"
include: "rules/bwa.smk"
include: "rules/bwa_mem2.smk"

rule all:
    input:
        compile_output_list(),
