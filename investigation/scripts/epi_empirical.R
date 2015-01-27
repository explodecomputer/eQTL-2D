# author: joseph.powell@uq.edu.au
# epi_empirical.R

# analysis to determine the empirical threshold

# Read in data 
# sig data
load('~/repo/eQTL-2D/investigation/data/sig_501_data.RData')	

# Genotype data
load('~/repo/eQTL-2D/data/geno.RData')

# Probe (residuals) data
load('~/repo/eQTL-2D/data/bsgs_egcut_fehr_data.RData')
bsgs <- phenlist[[1]]
rm(phenlist)
