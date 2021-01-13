library(tidyverse)
library(data.table)

args <- commandArgs(T)

discfile <- args[1]
replfile <- args[2]
varexp <- as.numeric(args[3])
i <- args[4]
sim <- as.numeric(args[5])
cis <- args[6]
sentinel <- args[7]
outfile <- args[8]
mapfile <- args[9]
famfile <- args[10]
signrepout <- args[11]

disc <- fread(discfile)  %>% as_tibble() %>% filter(V3 == 8)
print(nrow(disc))
disc$p4 <- pf(disc$V6, 4, disc$V4, low=FALSE)
disc$fdr <- p.adjust(disc$p4, "fdr")

repl <- fread(replfile) %>% as_tibble() %>% filter(V3 == 8)
print(nrow(repl))
repl$p4 <- pf(repl$V6, 4, repl$V4, low=FALSE)
repl$fdr <- p.adjust(repl$p4, "fdr")

map <- fread(mapfile) %>% as_tibble()

disc$rsid <- map$V2[disc$V2]
repl$rsid <- map$V2[repl$V2]
disc$chr <- map$V1[disc$V2]
repl$chr <- map$V1[repl$V2]
disc$pos <- map$V4[disc$V2]
repl$pos <- map$V4[repl$V2]


mer <- merge(disc, repl, by="V2")
print(head(mer))
print(dim(mer))

# Nsig significant in disc
# Lambda in disc
# Fstat corr in disc and repl
# nrep @ p<0.05
# nrep @ fdr<0.05

estlambda <- function(x)
{
	da <- qnorm(1 - x/2)
	# da <- qchisq(x, 1, low=FALSE)
	median(da, na.rm=TRUE)^2/qchisq(0.5, 1)
} 

estlambda2 <- function(F)
{
	median(F, na.rm=TRUE)/qchisq(0.5, 4) * 4
} 

res <- list()

res$cis <- cis
res$sentinel <- sentinel
res$varexp <- varexp
res$i <- i
res$sim <- sim

res$lambda_disc <- estlambda(disc$p4)
res$lambda_repl <- estlambda(repl$p4)
res$lambda_disc2 <- estlambda2(disc$V6)
res$lambda_repl2 <- estlambda2(repl$V6)

disc_bonf <- subset(disc, V3 == 8 & p4 < (0.05/nrow(disc))) %>%
arrange(p4) %>% 
filter(!duplicated(chr))
res$nsig_disc_bonf <- nrow(disc_bonf)

disc_bonf2 <- subset(disc, V3 == 8 & p4 < 1e-16) %>%
arrange(p4) %>% 
filter(!duplicated(chr))
res$nsig_disc_bonf2 <- nrow(disc_bonf2)

disc_fdr <- subset(disc, V3 == 8 & fdr < 0.05) %>%
arrange(fdr) %>% 
filter(!duplicated(chr))
res$nsig_disc_fdr <- nrow(disc_fdr)

repl_bonf <- subset(repl, V3 == 8 & p4 < (0.05/nrow(repl))) %>%
arrange(p4) %>% 
filter(!duplicated(chr))
res$nsig_repl_bonf <- nrow(repl_bonf)

repl_bonf2 <- subset(repl, V3 == 8 & p4 < 1e-16) %>%
arrange(p4) %>% 
filter(!duplicated(chr))
res$nsig_repl_bonf <- nrow(repl_bonf)

repl_fdr <- subset(repl, V3 == 8 & fdr < 0.05) %>%
arrange(fdr) %>% 
filter(!duplicated(chr))
res$nsig_repl_fdr <- nrow(repl_fdr)


cond_bonf <- subset(mer, V2 %in% disc_bonf$V2)
cond_bonf$fdr <- p.adjust(cond_bonf$p4.y, "fdr")

res$cond_bonf_bonf <- sum(cond_bonf$p4.y < (0.05 / nrow(cond_bonf)), na.rm=TRUE)
res$cond_bonf_bonf2 <- sum(cond_bonf$p4.y < (0.05 / 501), na.rm=TRUE)
res$cond_bonf_fdr <- sum(cond_bonf$fdr < 0.05, na.rm=TRUE)

cond_bonf2 <- subset(mer, V2 %in% disc_bonf2$V2)
cond_bonf2$fdr <- p.adjust(cond_bonf2$p4.y, "fdr")

res$cond_bonf2_bonf <- sum(cond_bonf2$p4.y < (0.05 / nrow(cond_bonf2)), na.rm=TRUE)
res$cond_bonf2_bonf2 <- sum(cond_bonf2$p4.y < (0.05 / 501), na.rm=TRUE)
res$cond_bonf2_fdr <- sum(cond_bonf2$fdr < 0.05, na.rm=TRUE)


cond_fdr <- subset(mer, V2 %in% disc_fdr$V2)
cond_fdr$fdr <- p.adjust(cond_fdr$p4.y, "fdr")

res$cond_fdr_bonf <- sum(cond_fdr$p4.y < (0.05 / nrow(cond_fdr)), na.rm=TRUE)
res$cond_fdr_bonf2 <- sum(cond_fdr$p4.y < (0.05 / 501), na.rm=TRUE)
res$cond_fdr_fdr <- sum(cond_fdr$fdr < 0.05, na.rm=TRUE)

res$fstat_cor <- cor(mer$V6.x, mer$V6.y, use="pair")

print(res)
save(res, file=outfile)


if(is.null(famfile))
{
	q()
}

if(nrow(disc_bonf) == 0)
{
	q()
}

# Do sign replication analysis if there are some significant hits

# extract data from discovery and replication

# perform noia on discovery and replication
# save decompositions for discovery and replication datasets


