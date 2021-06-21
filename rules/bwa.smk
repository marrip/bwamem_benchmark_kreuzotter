rule bwa:
    input:
        unpack(get_fastq),
        ref=config["reference"]["fasta"],
    output:
        temp("analysis_output/{unit}/bwa/{sample}_{threads}.bam"),
    log:
        "analysis_output/{unit}/bwa/{sample}_{threads}.log",
    benchmark:
        "analysis_output/{unit}/bwa/{sample}_{threads}.tsv",
    params:
        K=10000000,
        R="'@RG\\tID:{sample}_rg1\\tLB:lib1\\tPL:bar\\tSM:{sample}\\tPU:{sample}_rg1'",
    container:
        config["tools"]["common"]
    threads: get_threads
    message:
        "{rule}: Align {wildcards.sample} with {wildcards.threads} threads and sort records"
    shell:
        "(bwa mem "
        "-t {threads} "
        "-K {params.K} "
        "-R {params.R} "
        "{input.ref} "
        "{input.fwd} "
        "{input.rev} | "
        "samtools sort "
        "-@ {threads} "
        "-o {output} -) &> {log}"


rule index_bwa:
    input:
        "analysis_output/{unit}/bwa/{sample}_{threads}.bam",
    output:
        temp("analysis_output/{unit}/bwa/{sample}_{threads}.bai"),
    log:
        "analysis_output/{unit}/bwa/{sample}_{threads}_index_bam.log",
    container:
        config["tools"]["common"]
    message:
        "{rule}: Index {wildcards.sample} bam file"
    shell:
        "samtools index "
        "-b {input} "
        "{output} &> {log}"
