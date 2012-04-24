# Gib Hemani
# Create phenotype file for gcta

objfile <- commandArgs(T)[1] # famdat, phendat
id <- as.numeric(commandArgs(T)[2])
outfile <- commandArgs(T)[3]

load(objfile)

phen <- phendat[, id]

ids <- data.frame(famdat[,1], famdat[,2], phen)
write.table(ids, file=outfile, row=F, col=F, qu=F)

