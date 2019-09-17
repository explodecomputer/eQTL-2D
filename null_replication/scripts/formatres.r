library(tidyverse)
args <- commandArgs(T)

discfile <- args[1]
replfile <- args[2]
varexp <- as.numeric(args[3])
i <- as.numeric(args[4])
sim <- as.numeric(args[5])
cis <- args[6]
sentinel <- args[7]
outfile <- args[8]
mapfile <- args[9]


disc <- read_delim(discfile, "\t", col_names=FALSE)
disc$p4 <- pf(disc$X6, 4, disc$X4, low=FALSE)
disc$fdr <- p.adjust(disc$p4, "fdr")

repl <- read_delim(replfile, "\t", col_names=FALSE)
repl$p4 <- pf(repl$X6, 4, repl$X4, low=FALSE)
repl$fdr <- p.adjust(repl$p4, "fdr")

map <- read_delim(mapfile, "\t", col_names=FALSE)

disc$rsid <- map$X2[disc$X2]
repl$rsid <- map$X2[repl$X2]
disc$chr <- map$X1[disc$X2]
repl$chr <- map$X1[repl$X2]
disc$pos <- map$X4[disc$X2]
repl$pos <- map$X4[repl$X2]


mer <- merge(disc, repl, by="X2")

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
res$sim <- sim

res$lambda_disc <- estlambda(disc$p4)
res$lambda_repl <- estlambda(repl$p4)

disc_bonf <- subset(disc, X3 == 8 & p4 < (0.05/nrow(disc))) %>%
arrange(p4) %>% 
filter(!duplicated(chr))
res$nsig_disc_bonf <- nrow(disc_bonf)

disc_fdr <- subset(disc, X3 == 8 & fdr < 0.05) %>%
arrange(fdr) %>% 
filter(!duplicated(chr))
res$nsig_disc_fdr <- nrow(disc_fdr)

repl_bonf <- subset(repl, X3 == 8 & p4 < (0.05/nrow(repl))) %>%
arrange(p4) %>% 
filter(!duplicated(chr))
res$nsig_repl_bonf <- nrow(repl_bonf)

repl_fdr <- subset(repl, X3 == 8 & fdr < 0.05) %>%
arrange(fdr) %>% 
filter(!duplicated(chr))
res$nsig_repl_fdr <- nrow(repl_fdr)


cond_bonf <- subset(mer, X2 %in% disc_bonf$X2)
cond_bonf$fdr <- p.adjust(cond_bonf$p4.y, "fdr")

res$cond_bonf_bonf <- sum(cond_bonf$p4.y < 0.05 / nrow(cond_bonf), na.rm=TRUE)
res$cond_bonf_fdr <- sum(cond_bonf$fdr < 0.05, na.rm=TRUE)

cond_fdr <- subset(mer, X2 %in% disc_fdr$X2)
cond_fdr$fdr <- p.adjust(cond_fdr$p4.y, "fdr")

res$cond_fdr_bonf <- sum(cond_fdr$p4.y < 0.05 / nrow(cond_fdr), na.rm=TRUE)
res$cond_fdr_fdr <- sum(cond_fdr$fdr < 0.05, na.rm=TRUE)

res$fstat_cor <- cor(mer$X6.x, mer$X6.y, use="pair")

print(res)
save(res, file=outfile)
