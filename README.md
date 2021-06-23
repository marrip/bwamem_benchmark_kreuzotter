# bwamem_benchmark_kreuzotter

Benchmarking bwa and bwa-mem2 on slurm cluster

![Snakefmt](https://github.com/marrip/wgs_std_viper/actions/workflows/main.yaml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## :speech_balloon: Introduction

This workflow is designed to run **bwa** and **bwa-mem2** on a slurm cluster using
different thread counts to determine which algorithm performs better and what the
ideal thread count is. The results are plotted and the output `.bam` files are
compared to identify potential differences.

## :heavy_exclamation_mark: Dependencies

To run this workflow, the following tools need to be available:

![python](https://img.shields.io/badge/python-3.8-blue)

[![snakemake](https://img.shields.io/badge/snakemake-5.32.0-blue)](https://snakemake.readthedocs.io/en/stable/)

[![singularity](https://img.shields.io/badge/singularity-3.7-blue)](https://sylabs.io/docs/)

## :school_satchel: Preparations

### Sample data

1. Add all sample ids to `samples.tsv` in the column `sample`.
2. Add all sample data information to `units.tsv`. Each row represents a `fastq` file pair with
corresponding forward and reverse reads. Also indicate the sample id, run id and lane number.

### Reference data

1. You need a reference `.fasta` file to map your reads to. For the different tools to work, you also
need to prepare index files.

- Most of the required files for the human reference genome GRCh38 can be downloaded from
[google cloud](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0).
The download can be manually done using the browser or using `gsutil` via the command line:

```bash
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta /path/to/download/dir/
```

- If those resources are not available for your reference you may generate them yourself:

```bash
bwa index /path/to/reference.fasta
bwa-mem2 index /path/to/reference.fasta
```

2. Add the paths of the different files to the `config.yaml`. Index files should be
in the same directory as the reference `.fasta`.
3. Make sure that adapter sequences and docker container versions are correct.

## :white_check_mark: Testing

The workflow repository contains a small test dataset `.tests/integration` which can be run like so:

```bash
cd .tests/integration
snakemake -s ../../Snakefile -j1 --use-singularity
```

## :rocket: Usage

It is recommended to use a cluster profile and run something like:

```bash
snakemake -s /path/to/Snakefile --profile my-awesome-profile
```
