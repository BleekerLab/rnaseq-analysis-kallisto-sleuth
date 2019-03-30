"""
A Snakemake pipeline to go from mRNA-Seq reads to normalised transcript abundance estimates and differential expression
"""

############################
## Minimal Snakemake version
############################
#min_version("5.2.0")

### get liabraries
import pandas as pd

#############################################
## Configuration (parameters, samples, units)
#############################################

configfile: "config.yaml"

WORKING_DIR = config["workdir"]
RESULT_DIR  = config["resultdir"]

# get list of samples
units   = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

# directory that contains original fastq files
FQ_DIR  = config["fqdir"]
SAMPLES = units.index.get_level_values('sample').unique().tolist()
print(SAMPLES)

# Threads
THREADS = 10

##########################################
## Functions required to fetch input files
##########################################

def reads_are_SE(sample):
    """This function detect missing value in the column 2 of the units.tsv"""
    return pd.isnull(units.loc[(sample), "fq2"])

def get_fastq(wildcards):
	""" This function checks if sample is paired end or single end
	and returns a pair or single fastq file """
    if reads_are_SE(wildcards.sample):
        return units.loc[(wildcards.sample), ["fq1"]].dropna()
    else:
        return units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_trimmed(wildcards):
	""" This function checks if sample is paired end or single end
	and returns a pair or single filename """
    if reads_are_SE(wildcards.sample):
        return wildcards.sample + "_R1_trimmed.fq"
    else:
        return [wildcards.sample + "_R1_trimmed.fq", wildcards.sample + "_R2_trimmed.fq"]

##################
## Desired outputs
##################

KALLISTO = expand("results/kallisto/{samples}/abundance.tsv",samples=SAMPLES)
MASTERS = ["results/Snakefile","results/config.yaml","environment.yaml"]

rule all:
	input:
		KALLISTO,
		MASTERS
	message:"all done"

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

##############################################################################
## Kallisto (pseudo-alignment) analysis for transcriptome and custom databases
##############################################################################

rule estimate_transcript_abundance_using_kallisto:
    input:
        index = "index/kallisto_index.kidx",
        reads = get_trimmed
    output:
        "results/kallisto/{sample}/abundance.tsv"
    message:"computing {wildcards.sample} abundances using kallisto"
    params:
        sampleName      = "{sample}",
        outDir          = "results/kallisto/{sample}/",
        fragmentLength  = str(config["kallisto"]["fragment-length"]),
        sd              = str(config["kallisto"]["sd"]),
        bootstrap       = str(config["kallisto"]["bootstrap"])
    run:
        if reads_are_SE(params.sampleName):
            shell("mkdir -p results/kallisto/; \
            kallisto quant -i {input.index} -o {params.outDir} \
            --single -l {params.fragmentLength} -s {params.sd} \
            -b {params.bootstrap} \
            --threads {THREADS} \
            {input.reads}")
        else:
            shell("mkdir -p results/kallisto/; \
            kallisto quant -i {input.index} -o {params.outDir} \
            -b {params.bootstrap} \
            --threads {THREADS} \
            {input.reads}")


rule create_kallisto_index:
    input:
         fasta = config["kallisto"]["fasta_ref"]
    output:
         "index/kallisto_index.kidx"
    params:
        "kallisto_index.kidx"
    message:"creating kallisto index"
    shell:
        "kallisto index --make-unique -i {params} {input};"
        "mv {params} index/"

########################################################
## Read trimming of adapters, quality and max/min length
########################################################

rule trimmomatic:
    input:
        get_fastq
    output:
        fq1 = "{sample}_R1_trimmed.fq",
        fq2 = "{sample}_R2_trimmed.fq"
    message: "Trimming single-end {wildcards.sample} reads"
    log:
        RESULT_DIR + "logs/trimmomatic_se/{sample}.log"
    params :
        sampleName =                "{sample}",
        fq3 =                       "{sample}_R1_trimmed_unpaired.fq",
        fq4 =                       "{sample}_R2_trimmed_unpaired.fq",
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
        if reads_are_SE(params.sampleName):
            shell("trimmomatic SE {params.phred} -threads {THREADS} \
			{input} {output.fq1} \
			ILLUMINACLIP:{params.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} \
			LEADING:{params.LeadMinTrimQual} \
			TRAILING:{params.TrailMinTrimQual} \
			SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} \
			MINLEN:{params.minReadLen} CROP:{params.maxLen} ;\
			touch {output.fq2}")
        else:
            shell("trimmomatic PE {params.phred} -threads {THREADS} \
			{input} {output.fq1} {params.fq3} {output.fq2} {params.fq4} \
			ILLUMINACLIP:{params.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} \
			LEADING:{params.LeadMinTrimQual} \
			TRAILING:{params.TrailMinTrimQual} \
			SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} \
			MINLEN:{params.minReadLen} CROP:{params.maxLen}")
