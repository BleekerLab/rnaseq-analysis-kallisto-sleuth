#!/usr/bin/env python3

"""
This script is meant to subsample a single-end fastq file (zipped or unzipped) and randomly extract a specified number of reads. 
This is useful for tests etc. 
"""


###########
# Libraries
###########
import subprocess
import argparse # Read arguments passed by the command-line using the argparse library

######################
# Activate environment
######################

# activate the ad hoc conda environment for seqtk and argparse versioning
subprocess.call("source activate rnaseq-kallisto-sleuth",shell=True)


################################################
# Parse arguments given through the command-line
################################################
parser = argparse.ArgumentParser(description='Arguments to be passed to the seqtk command-line program.')

parser.add_argument("--fastq","-fq",type=str,help="The sequencing read fastq file to subsample",required=True)
parser.add_argument("--seed","-s",type=int,default=100,help="The seed used to randomly subsample")
parser.add_argument("--number","-n",type=int,default=100000,help="The number of reads to subsample (default = 100,000)")
parser.add_argument("--outfile","-o",type=str,default="subsampled.fastq",help="The name and path of the output subsampled file")


args = parser.parse_args()

###################################
# Have correct arguments been used?
###################################

# if fastq is gunzipped please unzip
# input can be gzipped or gunzipped, output is by defauld gunzipped

#############
# Subsampling
#############
seqtk_command = f"seqtk sample -s{args.seed} {args.fastq} {args.number} > {args.outfile}" 
subprocess.call(seqtk_command,shell=True)
