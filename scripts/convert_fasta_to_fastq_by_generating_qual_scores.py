"""
Convert FASTA to FASTQ file with a static
Usage:
$ ./fasta_to_fastq -i NAME.fasta -o NAME.fastq -m 30

Source: https://www.biostars.org/p/99886/

For more on Phred scores see https://www.drive5.com/usearch/manual/quality_score.html
"""

import argparse
from Bio import SeqIO
import random

# parse the command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i","--fasta_file", help = "The path to the input fasta file to convert")
parser.add_argument("-o","--fastq_file", help = "The path to the output fastq file file to write")
parser.add_argument("-m", "--min_quality", type = int, help = "The minimal Phred33 score quality value", default = 30)
parser.add_argument("-M", "--max_quality", type=int, help = "The maximal Phred33 score quality value ", default = 40)
args = parser.parse_args()


# Get inputs
fa_path = args.fasta_file
fq_path = args.fastq_file

# min and max values for Phred quality scores (Phred33)
min_qual = args.min_quality
max_qual = args.max_quality

# make fastq
with open(fa_path, "r") as fasta, open(fq_path, "w") as fastq:
    for record in SeqIO.parse(fasta, "fasta"):
        phred_score_qualities = random.choices(
        	list(range(min_qual, max_qual + 1, 1)),
            k = len(record)
            )
        record.letter_annotations["phred_quality"] = phred_score_qualities
        SeqIO.write(record, fastq, "fastq")