library(reshape2)
library(plyr)
library(xtable)


#=================================================================================================#
#=================================================================================================#


load("~/repo/eQTL-2D/analysis/interaction_list_meta_analysis.RData")
sig <- subset(meta, filter !=3)


#=================================================================================================#
#=================================================================================================#


# Results table 1

bsig <- subset(sig, pnest_meta > -log10(0.05/434))
bsig <- bsig[order(bsig$probegene), ]
dim(bsig)

tab <- subset(bsig, select=c(probegene, probechr, snp1, chr1, snp2, chr2, pnest, pnest_fehr, pnest_egcut, pnest_meta))
rownames(tab) <- 1:nrow(tab)

tab$snp1 <- paste(tab$snp1, " (", tab$chr1, ") ", sep="")
tab$snp2 <- paste(tab$snp2, " (", tab$chr2, ") ", sep="")
tab$probegene <- paste(tab$probegene, " (", tab$probechr, ") ", sep="")
tab <- subset(tab, select=-c(chr1, chr2, probechr))
tab

names(tab) <- c("Gene (chr.)", "SNP 1 (chr.)", "SNP 2 (chr.)", "Discovery", "Fehrmann", "EGCUT", "Combined replication")
xtable(tab, digits = c(0, 0, 0, 0, 2, 2, 2, 2))


#=================================================================================================#
#=================================================================================================#


# Supplementary tables


sigm <- sig
load("~/repo/eQTL-2D/analysis/interaction_list_replication_summary.RData")
dim(sig)
dim(sigm)

temp <- subset(sigm, select=c(code, pnest_meta))
sig <- merge(sig, temp, by="code", all.x=T)

a <- subset(sig, select=c(probegene, probename, probechr, snp1, chr1, position1, marginal_gene1, snp2, chr2, position2, marginal_gene2, pnest, pnest_fehr, pnest_egcut, pnest_meta))
dim(a)
head(a)

# bp diff for cis

a$dif <- NA
index <- a$chr1 == a$probechr & a$chr2 == a$probechr
a$dif[index] <- abs(a$position1 - a$position2)[index] / 1000000

a <- a[order(a$probegene), ]
head(a)

print(xtable(a, digits=c(rep(0, 12), 2, 2, 2, 2, 3)), include.rownames=FALSE)
write.csv(a, file="~/repo/eQTL-2D/analysis/results.csv", qu=F)
