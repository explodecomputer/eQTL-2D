library(tidyverse)
args <- commandArgs(T)

discfile <- args[1]
replfile <- args[2]
varexp <- as.numeric(args[3])
i <- as.numeric(args[4])
cis <- args[5]
sentinel <- args[6]
outfile <- args[7]
mapfile <- args[8]


disc <- read_delim(discfile, "\t", col_names=FALSE)
disc$p4 <- pf(disc$X6, 4, disc$X4, low=FALSE)

repl <- read_delim(replfile, "\t", col_names=FALSE)
repl$p4 <- pf(repl$X6, 4, repl$X4, low=FALSE)

mer <- merge(disc, repl, by="X2")

map <- read_delim(mapfile, "\t", col_names=FALSE)




# Nsig significant in disc
# Lambda in disc
# Fstat corr in disc and repl
# nrep @ p<0.05
# nrep @ fdr<0.05

estlambda <- function(x)
{
	da <- qchisq(x, 1, low=FALSE)
	median(da, na.rm=TRUE)/qchisq(0.5, 1)
} 

res <- list()

res$cis <- cis
res$sentinel <- sentinel
res$varexp <- varexp
res$i <- i

res$lambda_disc <- estlambda(disc$p4)
res$lambda_repl <- estlambda(repl$p4)

res$nsig_disc_p005 <- sum(disc$p4 < 0.05, na.rm=TRUE)
res$nsig_disc_fdr005 <- sum(p.adjust(disc$p4, "fdr") < 0.05, na.rm=TRUE)
res$nsig_disc_bonf <- sum(disc$p4 < (0.05/nrow(disc)), na.rm=TRUE)

res$nsig_repl_p005 <- sum(repl$p4 < 0.05, na.rm=TRUE)
res$nsig_repl_fdr005 <- sum(p.adjust(repl$p4, "fdr") < 0.05, na.rm=TRUE)
res$nsig_repl_bonf <- sum(repl$p4 < (0.05/nrow(repl)), na.rm=TRUE)

temp <- subset(mer, p4.x < (0.05/nrow(disc)))
temp$fdr <- p.adjust(temp$p4.y, "fdr")

res$cond1_p005 <- sum(temp$p4.y < 0.05, na.rm=TRUE)
res$cond1_fdr005 <- sum(temp$fdr < 0.05, na.rm=TRUE)

temp <- subset(mer, p.adjust(p4.x, "fdr") < 0.05)
temp$fdr <- p.adjust(temp$p4.y, "fdr")

res$cond2_p005 <- sum(temp$p4.y < 0.05, na.rm=TRUE)
res$cond2_fdr005 <- sum(temp$fdr < 0.05, na.rm=TRUE)

res$fstat_cor <- cor(mer$X6.x, mer$X6.y, use="pair")

print(res)
save(res, file=outfile)
