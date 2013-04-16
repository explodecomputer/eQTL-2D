# Gib Hemani
# Create .fam file for epiGPU

residfile <- commandArgs(T)[1]

residdat <- read.table(residfile, header=F)

phen <- residdat[, 6]

# normalise phenotype
# rank transform:

library(GenABEL)
phen <- rntransform(phen)

write.table(phen, file=paste(residfile, "phen", sep="."), row=F, col=F, qu=F)

