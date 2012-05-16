# Gib Hemani
# Fit mixed model to phenotype and return residuals

library(GenABEL)

datafile <- commandArgs(T)[1]
probe <- as.numeric(commandArgs(T)[2])
out <- commandArgs(T)[3]
load(datafile)
ex <- data.frame(phen = phendat[, probe+2], fac1 = covdat$V3, fac2 = covdat$V4)
rownames(ex) <- phendat$IID
h <- polygenic(phen ~ as.factor(fac1) + as.factor(fac2), kin=Amat, ex)
write.table(h$esth2, file=paste(out, "h2", sep="."), row=F, col=F, qu=F)
phenres <- rntransform(h$pgresidualY)
write.table(phenres, file=out, row.names=F, col.names=F, quote=F)

