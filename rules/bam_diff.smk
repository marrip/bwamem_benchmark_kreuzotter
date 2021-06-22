rule bam_diff:
    input:
        bam1="analysis_output/{unit}/bwa/{sample}_{threads}.bam",
        bam2="analysis_output/{unit}/bwa-mem2/{sample}_{threads}.bam",
    output:
        "analysis_output/{unit}/bam_diff/{sample}_{threads}.diff",
    log:
        "analysis_output/{unit}/bam_diff/{sample}_{threads}.log",
    container:
        config["tools"]["bamutil"]
    message:
        "{rule}: Compare {input.bam1} and {input.bam1}"
    shell:
        "bam diff "
        "--in1 {input.bam1} "
        "--in2 {input.bam2} "
        "--out {output} "
        "--noCigar "
        "--isize "
        "--flag "
        "--mate "
        "--mapQual &> {log} && "
        "touch {output}"
