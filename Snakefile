"""
A Snakemake pipeline to go from mRNA-Seq reads to normalised transcript abundance estimates and differential expression
"""

############################
## Minimal Snakemake version
############################
#min_version("5.2.0")

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
MASTERS = ["results/Snakefile","results/config.yaml","environment.yaml"]

rule all:
	input:
		KALLISTO,
		MASTERS
	message:"all done"
    #shell:
#        "rm -r {WORKING_DIR}"

################################
## Copy master files to results
##############################
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
    conda:
        "envs/kallisto.yaml"
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
        "kallisto_index.kidx"
    message:"creating kallisto index"
    conda:
        "envs/kallisto.yaml"
    shell:
        "kallisto index --make-unique -i {params} {input};"
        "mv {params} index/"

################
## Read trimming
################

rule trimmomatic:
    input:
        get_fastq
    output:
        fq1 = WORKING_DIR + "{sample}_R1_trimmed.fq.gz",
        fq2 = WORKING_DIR + "{sample}_R2_trimmed.fq.gz"
    message: "Trimming {wildcards.sample} reads"
    conda:
        "envs/trimmomatic.yaml"
    log:
        RESULT_DIR + "logs/trimmomatic_se/{sample}.log"
    params :
        sampleName =                "{sample}",
        fq1_unpaired =              WORKING_DIR + "{sample}_R1_trimmed_unpaired.fq",
        fq2_unpaired =              WORKING_DIR + "{sample}_R2_trimmed_unpaired.fq",
        seedMisMatches =            str(config['trimmomatic']['seedMisMatches']),
        palindromeClipTreshold =    str(config['trimmomatic']['palindromeClipTreshold']),
        simpleClipThreshhold =      str(config['trimmomatic']['simpleClipThreshold']),
        LeadMinTrimQual =           str(config['trimmomatic']['LeadMinTrimQual']),
        TrailMinTrimQual =          str(config['trimmomatic']['TrailMinTrimQual']),
        windowSize =                str(config['trimmomatic']['windowSize']),
        avgMinQual =                str(config['trimmomatic']['avgMinQual']),
        minReadLen =                str(config['trimmomatic']['minReadLength']),
        phred = 		            str(config["trimmomatic"]["phred"]),
        adapters =                  str(config["trimmomatic"]["adapters"]),
        maxLen =                    str(config["trimmomatic"]["maxLen"])
    run:
        if sample_is_single_end(params.sampleName):
            shell("trimmomatic SE {params.phred} -threads {threads} \
			{input} {output.fq1} \
			ILLUMINACLIP:{params.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} \
			LEADING:{params.LeadMinTrimQual} \
			TRAILING:{params.TrailMinTrimQual} \
			SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} \
			MINLEN:{params.minReadLen} CROP:{params.maxLen} 2>{log};\
			touch {output.fq2}")
        else:
            shell("trimmomatic PE {params.phred} -threads {threads} \
			{input} {output.fq1} {params.fq1_unpaired} {output.fq2} {params.fq2_unpaired} \
			ILLUMINACLIP:{params.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} \
			LEADING:{params.LeadMinTrimQual} \
			TRAILING:{params.TrailMinTrimQual} \
			SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} \
			MINLEN:{params.minReadLen} CROP:{params.maxLen} 2>{log}")
