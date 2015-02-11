# epi_empirical_analysis.R

# analysis and interpretation of the output from epi_empirical.R
# joseph.powell@uq.edu.au

# read in the functions
source("/Users/jpowell/repo/eQTL-2D/investigation/scripts/epi_empirical_analysis_fun.R")

# load data
load("/Users/jpowell/repo/eQTL-2D/investigation/data/investigation_data.RData")

# set the data dir
setwd("/Users/jpowell/repo/eQTL-2D/investigation/data/output/")

# list files
lf <- list.files()


##################################################################
##################################################################
##################################################################
# summary of the clean and total snps

genome_summary <- genome_sum.fun(lf)


##################################################################
##################################################################
##################################################################
# add the filter info and y/n in 30 sig

genome_summary <- filter_add.fun(sig, genome_summary)



##################################################################
##################################################################
##################################################################
# Determine the summary of the output

lambda <- summarize.fun(lf)



##################################################################
##################################################################
##################################################################
# Calculate the additive eQTL effects for each pair

add_test <- add_cal.fun(sig, bsgs, geno, bim)

