setwd("~/repo/eQTL-2D/filtering/marginal_lists")
a <- read.table("FDR_0_05_genotyped_eQTL.txt", header=T)
b <- read.table("dom.txt", he=T)

a <- subset(a, select=c(CHR, SNP, TRAIT, PVALUE))
b <- subset(b, select=c(CHR, SNP, V10, P))

names(a) <- names(b) <- c("chr", "snp", "probename", "pval")

a$eff <- "additive"
b$eff <- "dominant"

lapply(a, class)
lapply(b, class)
a$snp <- as.character(a$snp)
b$snp <- as.character(b$snp)
a$probename <- as.character(a$probename)
b$probename <- as.character(b$probename)

head(a)
head(b)

marginal_list <- rbind(a, b)
dim(marginal_list)
table(marginal_list$eff)

save(marginal_list, file="marginal_list.RData")
