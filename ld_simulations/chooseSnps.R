# R --no-save --args bsgs_target_snps.txt arichu_common_snps.txt.gz 100000


ar <- commandArgs(T)

chipsnpsfile <- ar[1]
snpfile <- ar[2]
nqtl <- as.numeric(ar[3])

chipsnps <- scan(chipsnpsfile, what="character")
allsnps <- scan(snpfile, what="character")

chooseSnps <- function(allsnps, chipsnps, nqtl)
{
	chip <- chipsnps[chipsnps %in% allsnps]
	left <- allsnps[! allsnps %in% chipsnps]
	qtls <- left[sample(1:length(left), nqtl, replace=FALSE)]
	l <- list()
	l$chip <- chip
	l$qtls <- qtls
	return(l)
}

snps <- chooseSnps(allsnps, chipsnps, nqtl)
lapply(snps, length)

gz1 <- gzfile("allsnps.txt.gz", "w")
write.table(c(snps$chip, snps$qtls), gz1, row=F, col=F, qu=F)
close(gz1)

gz2 <- gzfile("qtls.txt.gz", "w")
write.table(snps$qtls, gz2, row=F, col=F, qu=F)
close(gz2)
