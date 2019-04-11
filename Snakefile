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

configfile: "config.yaml"

WORKING_DIR = config["workdir"]
RESULT_DIR  = config["resultdir"]
FQ_DIR  = config["fqdir"]

########################
# Samples and conditions
########################
samples = pd.read_csv("samples.tsv", dtype=str,index_col=0,sep="\t")
SAMPLES = samples.index.get_level_values('sample').unique().tolist()

############################
## Input functions for rules
############################

def sample_is_single_end(sample):
    """This function detect missing value in the column 2 of the units.tsv"""
    if "fq2" not in samples.columns:
        return True
    else:
        return pd.isnull(samples.loc[(sample), "fq2"])

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
MASTERS  = ["results/Snakefile","results/config.yaml","environment.yaml"]
SLEUTH   = [RESULT_DIR + "sleuth_object.Rds",RESULT_DIR + "abundance_wide.tsv"]

rule all:
	input:
		KALLISTO,
		MASTERS,
		SLEUTH
	message:"all done"
    #shell:
#        "rm -r {WORKING_DIR}"


###############################
## Copy master files to results
###############################
rule copy_master_files_to_results:
    input:
        "Snakefile",
        "config.yaml",
        "environment.yaml",
        fasta = config["kallisto"]["fasta_ref"]
    output:
        RESULT_DIR + "Snakefile",
        RESULT_DIR + "config.yaml",
        RESULT_DIR + "environment.yaml"
    message:"copy master files to ./results/"
    shell:
        "cp {input} results/"

#######################
## create sleath object
#######################

rule run_sleuth:
    input:
        expand(RESULT_DIR + "kallisto/{samples}/abundance.tsv",samples=SAMPLES)
    output:
        sleuth_object     = RESULT_DIR + "sleuth_object.Rds",
        normalized_counts = RESULT_DIR + "abundance_wide.tsv"
    params:
        sample_file              = config["samples"],
        input_directory          = RESULT_DIR + "kallisto",
        p_value                  = config["sleuth"]["p_value"],
        output_directory         = RESULT_DIR
    threads: 10
    shell:
        "Rscript --vanilla scripts/sleuth_analysis.R -i {params.input_directory} "
        "-s {params.sample_file} "
        "-c {threads} "
        "-p {params.p_value} "
        "-o {params.output_directory} "


##############################################################################
## Kallisto (pseudo-alignment) analysis for transcriptome and custom databases
##############################################################################

rule estimate_transcript_abundance_using_kallisto:
    input:
        index = WORKING_DIR + "index/kallisto_index.kidx",
        reads = get_trimmed
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
        if sample_is_single_end(params.sampleName):
            shell("mkdir -p results/kallisto/; \
            kallisto quant -i {input.index} -o {params.outDir} \
            --single -l {params.fragmentLength} -s {params.sd} \
            -b {params.bootstrap} \
            --threads {threads} \
            {input.reads}")
        else:
            shell("mkdir -p results/kallisto/; \
            kallisto quant -i {input.index} -o {params.outDir} \
            -b {params.bootstrap} \
            --threads {threads} \
            {input.reads}")


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
        fq1 = WORKING_DIR + "{sample}_R1_trimmed.fq.gz",
        fq2 = WORKING_DIR + "{sample}_R2_trimmed.fq.gz",
        html = RESULT_DIR + "fastp/{sample}.html"
    message:"trimming {wildcards.sample} reads"
    threads: 10
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        sampleName = "{sample}",
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    run:
        if sample_is_single_end(params.sampleName):
            shell("fastp --thread {threads}  --html {output.html} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --in1 {input} --out1 {output} \
            2>{log}; \
			touch {output.fq2}")
        else:
            shell("fastp --thread {threads}  --html {output.html} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --detect_adapter_for_pe \
            --in1 {input[0]} --in2 {input[1]} --out1 {output.fq1} --out2 {output.fq2}; \
            2>{log}")
