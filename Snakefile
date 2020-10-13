    """
A Snakemake pipeline to go from mRNA-Seq reads to normalised transcript abundance estimates and differential expression
"""
from snakemake.utils import min_version

############################
## Minimal Snakemake version
############################
min_version("5.2.0")

############
## Libraries
############

import pandas as pd

###############
# Configuration
###############

configfile: "config/config.yaml"

WORKING_DIR = config["workdir"]
RESULT_DIR  = config["resultdir"]
FQ_DIR  = config["fqdir"]

########################
# Samples and conditions
########################
samples = pd.read_csv("config/samples.tsv", dtype=str, sep="\t").set_index("sample", drop=False)
SAMPLES = samples.index.tolist()

############################
## Input functions for rules
############################

def sample_is_single_end(sample):
    """This function detect missing value in the column 2 of the units.tsv"""
    return pd.isnull(samples.loc[sample,"fq2"])

def get_fastq(wildcards):
	""" This function checks if the sample has paired end or single end reads
	and returns 1 or 2 names of the fastq files """
	if sample_is_single_end(wildcards.sample):
		return samples.loc[(wildcards.sample), ["fq1"]].dropna()
	else:
		return samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_trimmed(wildcards):
	""" This function checks if sample is paired end or single end
	and returns 1 or 2 names of the trimmed fastq files """
	if sample_is_single_end(wildcards.sample):
		return WORKING_DIR + wildcards.sample + "_R1_trimmed.fq.gz"
	else:
		return [WORKING_DIR + wildcards.sample + "_R1_trimmed.fq.gz", WORKING_DIR + wildcards.sample + "_R2_trimmed.fq.gz"]

##################
## Desired outputs
##################

KALLISTO = expand(RESULT_DIR + "kallisto/{samples}/abundance.tsv",samples=SAMPLES)
TIDY_COUNTS = RESULT_DIR + "abundance_tidy.tsv"

rule all:
	input:
		KALLISTO,
		TIDY_COUNTS
	message:"all done"

###############################
## Rules
###############################

##########################################################
## Sleuth outputs: scaled counts and differential analysis
##########################################################

rule get_scaled_counts:
    input:
        expand(RESULT_DIR + "kallisto/{samples}/abundance.tsv",samples=SAMPLES)
    output:
        normalized_counts = RESULT_DIR + "abundance_tidy.tsv"
    params:
        sample_file              = config["samples"],
        input_directory          = RESULT_DIR + "kallisto",
        output_directory         = RESULT_DIR
    threads: 10
    shell:
        "Rscript --vanilla scripts/sleuth_analysis.R "
        "-i {params.input_directory} "
        "-s {params.sample_file} "
        "-c {threads} "
        "-o {params.output_directory} "


##############################################################################
## Kallisto (pseudo-alignment) analysis for transcriptome and custom databases
##############################################################################

rule estimate_transcript_abundance_using_kallisto:
    input:
        index = WORKING_DIR + "index/kallisto_index.kidx",
        fq1 = WORKING_DIR + "{sample}_R1_trimmed.fq.gz",
        fq2 = WORKING_DIR + "{sample}_R2_trimmed.fq.gz"
    output:
        RESULT_DIR + "kallisto/{sample}/abundance.tsv"
    message:"computing {wildcards.sample} abundances using kallisto"
    threads: 10
    params:
        sampleName      = "{sample}",
        outDir          = "results/kallisto/{sample}/",
        fragmentLength  = str(config["kallisto"]["fragment-length"]),
        sd              = str(config["kallisto"]["sd"]),
        bootstrap       = str(config["kallisto"]["bootstrap"])
    run:
        if sample_is_single_end(wildcards.sample):
            shell("mkdir -p results/kallisto/; \
            kallisto quant -i {input.index} -o {params.outDir} \
            --single -l {params.fragmentLength} -s {params.sd} \
            -b {params.bootstrap} \
            --threads {threads} \
            {input.fq1}")
        else:
            shell("mkdir -p results/kallisto/; \
            kallisto quant -i {input.index} -o {params.outDir} \
            -b {params.bootstrap} \
            --threads {threads} \
            {input.fq1} {input.fq2}")


rule create_kallisto_index:
    input:
        fasta = config["kallisto"]["fasta_ref"]
    output:
         WORKING_DIR + "index/kallisto_index.kidx"
    params:
        WORKING_DIR + "index/kallisto_index.kidx"
    message:"creating kallisto index"
    shell:
        "kallisto index --make-unique -i {params} {input};"

################
## Read trimming
################
rule fastp:
    input:
        get_fastq
    output:
        fq1 = temp(WORKING_DIR + "{sample}_R1_trimmed.fq.gz"),
        fq2 = temp(WORKING_DIR + "{sample}_R2_trimmed.fq.gz"),
        html = RESULT_DIR + "fastp/{sample}.html",
        json = RESULT_DIR + "fastp/{sample}.json",
    message:"trimming {wildcards.sample} reads"
    threads: 10
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        sample_name = "{sample}",
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    run:
        if sample_is_single_end(params.sample_name): # single end
            shell("fastp --thread {threads} --html {output.html} --json {output.json} --qualified_quality_phred {params.qualified_quality_phred}  --in1 {input} --out1 {output} 2>{log}")
            shell("touch {output.fq2}")
        else:  
            shell("fastp --thread {threads} --html {output.html} --json {output.json} --qualified_quality_phred {params.qualified_quality_phred}  --in1 {input[0]} --in2 {input[1]} --out1 {output.fq1} --out2 {output.fq2} 2>{log}")  




