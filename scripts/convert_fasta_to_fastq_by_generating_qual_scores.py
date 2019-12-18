"""
Convert FASTA to FASTQ file with a static
Usage:
$ ./fasta_to_fastq NAME.fasta NAME.fastq

Source: https://www.biostars.org/p/99886/

For more on Phred scores see https://www.drive5.com/usearch/manual/quality_score.html
"""

import sys, os
import argparse
from Bio import SeqIO

# Get inputs
fa_path = sys.argv[1]
fq_path = sys.argv[2]

# min and max values for Phred quality scores (Phred33)
min_qual = sys.argv[3]
max_qual = sys.argv[4]

# make fastq
with open(fa_path, "r") as fasta, open(fq_path, "w") as fastq:
    for record in SeqIO.parse(fasta, "fasta"):
    	quality = random.randrange(min_qual, max_qual, step = 1)
        record.letter_annotations["phred_quality"] = [quality] * len(record)
        SeqIO.write(record, fastq, "fastq")