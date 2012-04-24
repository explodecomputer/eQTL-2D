# Gib Hemani
# Create phenotype file for gcta

objfile <- commandArgs(T)[1] # famdat, phedat
id <- as.numeric(commandArgs(T)[2]

phen <- phedat[, id]

ids <- data.frame(famdat[,1], famdat[,2], phen)
write.table(ids, file=paste("gcta.phen", id, sep="."), row=F, col=F, qu=F)

