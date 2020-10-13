# A Snakemake RNA-Seq pipeline with Kallisto and Sleuth

A snakemake pipeline for the analysis of RNA-seq data that makes use of [Kallisto and sleuth](https://scilifelab.github.io/courses/rnaseq/labs/kallisto).

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Miniconda](https://img.shields.io/badge/miniconda-blue.svg)](https://conda.io/miniconda)

<!-- MarkdownTOC autolink="True" -->

- [Overview](#overview)
	- [Aim](#aim)
	- [Outputs](#outputs)
	- [Content](#content)
- [Installation](#installation)
	- [Download or clone the Github repository](#download-or-clone-the-github-repository)
	- [\(Option 1\) Installing and activating a virtual environment](#option-1-installing-and-activating-a-virtual-environment)
	- [\(Option 2\) Using a Docker image](#option-2-using-a-docker-image)
- [Configuration](#configuration)
	- [Configuration file \(config.yaml\)](#configuration-file-configyaml)
	- [Experimental design \(samples.tsv\)](#experimental-design-samplestsv)
	- [Test files](#test-files)
- [Snakemake execution](#snakemake-execution)
	- [\(Option 1\) Run within the conda environment](#option-1-run-within-the-conda-environment)
	- [\(Option 2\) Run within a Docker container](#option-2-run-within-a-docker-container)
	- [Cluster execution](#cluster-execution)
- [Graph of jobs](#graph-of-jobs)

<!-- /MarkdownTOC -->


# Overview 
## Aim
To perform the _pseudo-alignment_ steps of RNA-seq (Illumina) reads to a transcriptome reference and output individual Kallisto estimates and produce a file of the transcript scaled counts.  

## Outputs
This pipeline analyses the raw RNA-seq data and produces:
1. A file containing normalized counts.
2. Individual Kallisto estimates that can be used for differential expression with Sleuth. 

## Content
- `Snakefile`: a master file that contains the desired outputs and the rules to generate them from the input files.
- `config/config.yaml`: the configuration files making the Snakefile adaptable to any input files, transcriptome and parameters for the rules.
- `fastq/`: This folder should contain single-end or paired-end reads, or a mixture of paired and single-end reads in fastq format.
- `envs/`: a folder containing the environments needed for the conda package manager. If run with the `--use-conda` command, Snakemake will install the necessary softwares and packages using the conda environment files.
- `config/samples.tsv`:  a file containing information about the names, the paths and the conditions of the samples used.
**This file has to be adapted to your sample names before running the pipeline**.

# Installation

## Download or clone the Github repository
You will need a local copy of the `rnaseq-analysis-kallisto-sleuth` on your machine. You can either:
- use git in the shell: `git clone git@github.com:BleekerLab/rnaseq-analysis-kallisto-sleuth.git`
- click on "Clone or download" and select `download`

## (Option 1) Installing and activating a virtual environment
1. Download [Miniconda3](https://docs.conda.io/en/latest/miniconda.html)for your system (Windows, Mac OS X, Linux) that will act as the package/software manager.
2. Create a virtual environment named `kallisto` using the `environment.yaml` file with the following command: `conda env create --file environment.yaml`.
3. Then, activate this virtual environment with `source activate kallisto` or  `conda activate kallisto` (with `conda =>4.5.0`).

You should have everything you need.

## (Option 2) Using a Docker image

A custom-made Docker image is available [on DockerHub](https://hub.docker.com/r/bleekerlab/rnaseq-analysis-kallisto-sleuth?utm_source=docker4mac_2.3.0.4&utm_medium=repo_open&utm_campaign=referral). It contains all required softwares and packages (including Snakemake).

To use it:
1. Make sure Docker is available on your machine. [See instructions](https://docs.docker.com/get-docker/).
2. In your favorite Shell, pull it: `docker pull bleekerlab/rnaseq-analysis-kallisto-sleuth:4.7.12`

# Configuration

## Configuration file (config.yaml)
Make sure you have changed the parameters in the `config/config.yaml` file that specifies where to find:
- the sample data file `samples.tsv`
- the genomic and transcriptomic reference fasta files
- various parameters for certain softwares etc.    
This file is used so the `Snakefile` does not need to be changed when locations or parameters need to be changed.

## Experimental design (samples.tsv)
To get the right reads to the samples the `config/sample.tsv` needs to contain certain features.
- the column name of the column comtaining the sample names needs to be 'sample'
- the column names of of the columns containg the forward and reverse reads need to be 'fq1' and 'fq2'
- the column containing conditions, genotypes, treatment, etc is free of choise.

Here is an example of a file for an experiment containing paired-end reads:

| sample   | treatment | fq1 | fq2 |
| ------- | ---------- |-----|-----|
| sample1 | control | readsS1_R1.fastq | readsS1_R2.fastq |
| sample2 | control | readsS2_R1.fastq | readsS2_R2.fastq |
| sample3 | treated | readsS3_R1.fastq | readsS3_R2.fastq |
| sample4 | treated | readsS4_R1.fastq | readsS4_R2.fastq |


In case of an experiment containing only single-end reads, the column 'fq2' can be omitted.
The 'sample.tsv' will then look something like this:

| sample   | treatment | fq1 |
| ------- | ---------- |-----|
| sample1 | control | readsS1.fastq |
| sample2 | control | readsS2.fastq |
| sample3 | treated | readsS3.fastq |
| sample4 | treated | readsS4.fastq |


If the experiment contains both single and paired end reads, it should be something like this:

| sample   | condition | fq1 | fq2 |
| ------- | ---------- |-----|-----|
| sample1 | control | readsS1_R1.fastq | readsS1_R2.fastq |
| sample2 | control | readsS2.fastq | |
| sample3 | treated | readsS3_R1.fastq | readsS3_R2.fastq |
| sample4 | treated | readsS4.fastq | |

## Test files
A collection of test files are available on the [Zenodo archive here](https://doi.org/10.5281/zenodo.4085315). Create a `.test` folder and move the fastq files there. 

# Snakemake execution
The Snakemake pipeline/workflow management system reads a master file (often called `Snakefile`) to list the steps to be executed and defining their order. It has many rich features. Read more [here](https://snakemake.readthedocs.io/en/stable/).

## (Option 1) Run within the conda environment
1. Activate the `kallisto` conda environment: `conda activate kallisto`
2. Run with `snakemake --cores 1` or `snakemake -np` for a dry run.

## (Option 2) Run within a Docker container 
1. Place yourself in the `rnaseq-analysis-kallisto-sleuth/` folder.
2. Run the container by linking your current working directory within the `/home/snakemake/` folder inside of the container. `docker run --rm -v $PWD:/home/snakemake/ bleekerlab/rnaseq-analysis-kallisto-sleuth:4.7.12`. 

The docker run command triggers the "snakemake" command. You can add any Snakemake options after that:
* `docker run --rm -v $PWD:/home/snakemake/ bleekerlab/rnaseq-analysis-kallisto-sleuth:4.7.12 --cores N` where N is the number of cores.
* `docker run --rm -v $PWD:/home/snakemake/ bleekerlab/rnaseq-analysis-kallisto-sleuth:4.7.12 -np` for a dry run.
..etc...

## Cluster execution
For cluster execution, please refer to the [Snakemake reference](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution).

# Graph of jobs
![dag.png](./dag.png)
