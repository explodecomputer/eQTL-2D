# Read in pedigree information and phenotype data
# Correct phenotype for polygenic effects
# Output: .fam file with residuals from polygenic effects

library(pedigree)

objfile <- commandArgs(T)[1] # contains Ainv, fam, phedat
id <- as.numeric(commandArgs(T)[4]


phedat <- read.csv(phefile, header=F)

phen <- phedat[, id]

pheres <- Ainv %*% y

fam <- read.table(famfile, header=T)
fam[, 6] <- pheres
write.table(fam, file=paste(famfile, id, sep=""), row.names=F, col.names=F, quote=F)



