# Gib Hemani
# Create .fam file for epiGPU

famfile <- commandArgs(T)[1]
residfile <- commandArgs(T)[2]
id <- as.numeric(commandArgs(T)[3]

famdat <- read.table(fam, header=F)
residdat <- read.table(residfile, header=?)

famdat[,6] <- residdat[, ?]

write.table(famdat, file=paste(famfile, id, sep="."), row=F, col=F, qu=F)

