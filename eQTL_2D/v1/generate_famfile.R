# Read in pedigree information and phenotype data
# Correct phenotype for polygenic effects
# Output: .fam file with residuals from polygenic effects

library(pedigree)

objfile <- commandArgs(T)[1] # fam, phedat
id <- as.numeric(commandArgs(T)[4]


# make pheno file for gcta (include covariates)
# run gcta and read in residuals
# make .fam file with residuals
# exit
# script will run epiGPU

phen <- phedat[, id]

fam <- read.table(famfile, header=T)
fam[, 6] <- pheres
write.table(fam, file=paste(famfile, id, sep=""), row.names=F, col.names=F, quote=F)



