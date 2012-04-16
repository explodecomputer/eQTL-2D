# Read in pedigree information and phenotype data
# Correct phenotype for polygenic effects
# Output: .fam file with residuals from polygenic effects

phefile <- commandArgs(T)[1]
famfile <- commandArgs(T)[2]
id <- as.numeric(commandArgs(T)[3]

phedat <- read.table(phefile, header=F)

phen <- 


