suppressPackageStartupMessages(library("sleuth"))
suppressPackageStartupMessages(library("tidyverse"))

# capture the command-line arguments after --args (e.g. the shortstack results directory)
args <- commandArgs(trailingOnly = TRUE)     

# renaming for meaningful names
kallisto_directory = args[1] # indicates where are stored the kallisto results (one directory per sample)
number_of_cores = args[2]    # indicates how many cores we are using
design_file = args[3]        # design file that has both sample names and experimental conditions 
outdir = args[4]             # path to a directory where result tables will be saved
pval.signif = args[5]        # significance p-value threshold e.g. 0.05

samples = list.dirs(kallisto_directory,full.names = F,recursive = F)
kallisto_files = file.path(kallisto_directory,samples)

# read the experimental design file
# keep only samples and condition column
samples2condition = read.table(design_file,header = T,stringsAsFactors = F)
samples2condition = samples2condition[,1:2]

# add the path to kallisto result files
samples2condition = mutate(samples2condition,path = kallisto_files)
samples2condition

so <- sleuth_prep(samples2condition,num_cores = number_of_cores,extra_bootstrap_summary=FALSE)

###########################################################
# extract counts and perform differential expression tests
##########################################################
# outputs tidy format and human-readable wide format 
abundance.res.tidy = kallisto_table(so,use_filtered = TRUE,normalized = TRUE)
abundance.res.tidy = select(abundance.res.tidy,-tpm,-eff_len,-condition,-len)
abundance.res.wide = spread(data = abundance.res.tidy,key = target_id,value = est_counts)

# differential tests
#so <- sleuth_fit(so, ~condition, 'full')
#so <- sleuth_fit(so, ~1, 'reduced')
#o <- sleuth_lrt(so,"reduced","full")

#sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) %>%
 # dplyr::filter(qval <= pval.signif)

##############
# Save results
#############
# write tables 
write.table(x = abundance.res.tidy,file = file.path(outdir,"abundance_tidy.tsv"),quote = F,sep = "\t",row.names = F)
write.table(x = abundance.res.wide,file = file.path(outdir,"abundance_wide.tsv"),quote = F,sep = "\t",row.names = F)
#write.table(x = sleuth_table,file = file.path(outdir,"differential.tsv"),quote = F,sep = "\t",row.names = F)

# save R object
saveRDS(so,file = file.path(outdir,"sleuth_object.Rds"))