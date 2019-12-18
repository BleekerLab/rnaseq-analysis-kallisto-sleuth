library(polyester)
library(Biostrings)

# generate fold changes
fold_changes = matrix(data = c(4,4,rep(1,18),1,1,4,4,rep(1,16)), nrow=20)

# FASTA annotation
fasta_file = "../refs/TAIR10_cdna_20110103_representative_gene_model_updated.fasta"
fasta = readDNAStringSet(fasta_file)

# subset the FASTA file to first 20 transcripts
small_fasta = fasta[1:20]
writeXStringSet(small_fasta, '../refs/TAIR10_small.fa')

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
reads_per_transcript = round(width(small_fasta) / 100 * 20)

# simulation call:
simulate_experiment(fasta = '../refs/TAIR10_small.fa', 
                    reads_per_transcript = reads_per_transcript,         # here 20X coverage 
                    num_reps=c(4,4),                                   # number of biol replicates per group
                    fold_changes = fold_changes,                         # matrix of fold changes
                    size = NULL,                                        # negative binomial size (defaults to reads_per_transcript * fold_changes / 3)
                    paired = FALSE,
                    strand_specific = TRUE,
                    outdir='../data/simulated_reads/') 