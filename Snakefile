"""
A Snakemake pipeline to go from mRNA-Seq reads to normalised transcript abundance estimates and differential expression
"""

############################
## Minimal Snakemake version
############################
#min_version("5.2.0")

#############################################
## Configuration (parameters, samples, units)
#############################################

configfile: "config.yaml"

WORKING_DIR = config["workdir"]
RESULT_DIR = config["resultdir"]

# directory that contains original fastq files
FQ_DIR = config["fqdir"]
SAMPLES, = glob_wildcards(FQ_DIR + "{sample}.fq.gz")
print(SAMPLES)
# read length parameters
MIN_LEN = 25
MAX_LEN = 100

# Threads
THREADS = 10

####################
## Desired outputs
####################
KALLISTO = expand("results/kallisto/{samples}/abundance.tsv",samples=SAMPLES)
MASTERS = ["results/Snakefile","results/config.yaml","environment.yaml"]

rule all:
	input:
		KALLISTO,
		MASTERS
	message:"all done"

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
        index = "index/kallisto_index.kidx",
        reads = "trimmed/{sample}_qc.fq"
    output:
        "results/kallisto/{sample}/abundance.tsv"
    message:"computing {wildcards.sample} abundances using kallisto"
    params:
        outDir          = "results/kallisto/{sample}/",
        fragmentLength  = str(config["kallisto"]["fragment-length"]),
        sd              = str(config["kallisto"]["sd"]),
        bootstrap       = str(config["kallisto"]["bootstrap"])
    log:"results/kallisto/{sample}/log.txt"
    shell:
        "mkdir -p results/kallisto/;"
        "kallisto quant -i {input.index} -o {params.outDir} "
        "--single -l {params.fragmentLength} -s {params.sd} "			# average fragment length of 180 nts
        "-b {params.bootstrap} "					# number of bootstraps
        "--threads {THREADS} "
        "{input.reads} 2>{log}"

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

################
## Read trimming
################

rule trimmomatic_se:
    input:
        FQ_DIR + "{sample}.fq.gz"
    output:
        "trimmed/{sample}_qc.fq"
    message: "Trimming single-end {wildcards.sample} reads"
    log:
        RESULT_DIR + "logs/trimmomatic_se/{sample}.log"
    params :
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
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        "trimmomatic SE {params.phred} -threads {THREADS} "
        "{input} "
        "{output} "
        "ILLUMINACLIP:{params.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} "
        "LEADING:{params.LeadMinTrimQual} "
        "TRAILING:{params.TrailMinTrimQual} "
        "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
        "MINLEN:{params.minReadLen} "
        "CROP:{params.maxLen}"
