suppressPackageStartupMessages(library("sleuth"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("optparse"))

# command-line arguments to provide
option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default="results/kallisto", help="directory where the sample directories containing the abundance files are located", metavar="character"),
  make_option(c("-s", "--sample_file"), type="character", default="samples.tsv", help="sample files used to get conditions for DESEq2 model fit", metavar="character"),
  make_option(c("-c", "--number_of_cores"), type="integer", default=1, help="sample files used to get conditions for DESEq2 model fit", metavar="integer"),
  make_option(c("-p", "--p_vaue"), type="double", default=0.05, help="maximum p_value to be significant", metavar="double"),
  make_option(c("-o", "--outdir"), type="character", default="results", help="where to place differential expression files", metavar="character")
)

# parse the command-line arguments and pass them to a list called 'opt'
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# get list of names of sample directories containing the abandance.csv files
samples <- list.dirs(opt$input_dir, full.names = F, recursive = F)
kallisto_files <- file.path(opt$input_dir, samples)

# read the experimental design file
# keep only samples and condition column
samples2condition = read.table(opt$sample_file,header = T,stringsAsFactors = F, sep = "\t")
samples2condition = samples2condition[,1:2]

# add the path to kallisto result files
samples2condition = mutate(samples2condition,path = kallisto_files)

so <- sleuth_prep(samples2condition, num_cores = opt$number_of_cores, extra_bootstrap_summary=FALSE)

###########################################################
# extract counts and perform differential expression tests
##########################################################
# outputs tidy format and human-readable wide format
abundance.res.tidy = kallisto_table(so,use_filtered = TRUE,normalized = TRUE)
abundance.res.tidy = select(abundance.res.tidy,-tpm,-eff_len,-condition,-len)
abundance.res.wide = t(spread(data = abundance.res.tidy,key = target_id, value = est_counts))

# differential tests
# so <- sleuth_fit(so, ~condition, 'full')
# so <- sleuth_fit(so, ~1, 'reduced')
# so <- sleuth_lrt(so,"reduced","full")


#sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_significant <- dplyr::filter(sleuth_table, qval <= opt$p_vaue)

##############
# Save results
#############
# write tables
write.table(x = abundance.res.tidy, file = file.path(opt$outdir, "abundance_tidy.tsv"), quote = F,sep = "\t", row.names = F)
write.table(x = abundance.res.wide, file = file.path(opt$outdir, "abundance_wide.tsv"), quote = F,sep = "\t", row.names = T, col.names = F)
#write.table(x = sleuth_table,file = file.path(outdir,"differential.tsv"),quote = F,sep = "\t",row.names = F)

# save R object
saveRDS(so, file = file.path(opt$outdir,"sleuth_object.Rds"))
