library(tidyverse)
args <- commandArgs(T)

resfile <- "../data/discready.res"
bimfile <- "../data/discready.bim"

resfile <- args[1]
bimfile <- args[2]
outfile <- args[3]


res <- read_delim(resfile, "\t", col_names=FALSE)
res$p4 <- pf(res$X6, 4, res$X4, low=FALSE)
res$fdr4 <- p.adjust(res$p4, "fdr")
res <- subset(res, fdr4 < 0.5)

bim <- read_delim(bimfile, "\t", col_names=FALSE)

res$X2 <- res$X2 + 1
rsid <- bim$X2[res$X2]

write.table(rsid, file=outfile, row=F, col=F, qu=F)

