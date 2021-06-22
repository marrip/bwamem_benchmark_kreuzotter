import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.4.0")

### Set and validate config file


configfile: "config.yaml"


validate(config, schema="../schemas/config.schema.yaml")

### Read and validate samples file

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = pd.read_table(config["units"], dtype=str).set_index(
    ["sample", "unit"], drop=False
)
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    sample="|".join(samples.index),
    unit="|".join(units["unit"]),


### Functions


def get_fastq(wildcards):
    fastqs = units.loc[
        (wildcards.sample, wildcards.unit), ["fq1", "fq2"]
    ].dropna()
    return {"fwd": fastqs.fq1, "rev": fastqs.fq2}


def get_threads(wildcards):
  return config["threads"][wildcards.threads]


def get_bm_files():
    return expand(
        "analysis_output/{unit}/{tool}/{sample}_{threads}.tsv",
        unit = units["unit"],
        tool = ["bwa", "bwa-mem2"],
        sample = samples.index,
        threads = config["threads"],
        )


def compile_output_list():
    return expand(
        "analysis_output/{unit}/plot.pdf",
        unit = units["unit"],
        ) + expand(
        "analysis_output/{unit}/bam_diff/{sample}_{threads}.diff",
        unit = units["unit"],
        sample = samples.index,
        threads = config["threads"],
        )
