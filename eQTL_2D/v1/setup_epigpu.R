# Gib Hemani
# Create .fam file for epiGPU

famfile <- commandArgs(T)[1]
residfile <- commandArgs(T)[2]

famdat <- read.table(famfile, header=F)
residdat <- read.table(residfile, header=F)

famdat[,6] <- residdat[, 6]

write.table(famdat, file=paste(residfile, "fam", sep="."), row=F, col=F, qu=F)

