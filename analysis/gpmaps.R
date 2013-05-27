library(lattice)
library(latticeExtra)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(plyr)

plot3dHairballsRep <- function(sig, pg, cissnp, z=45)
{
	# make table of mean phenotypes
	# AA cis vs 1:nalleles trans
	# Aa cis vs 1:nalleles trans
	# aa cis vs 1:nalleles trans

	s <- subset(sig, probegene==pg)
	pn <- s$probename[1]
	snps <- unique(c(s$snp1, s$snp2))
	stopifnot(cissnp %in% snps)
	tsnps <- snps[! snps %in% cissnp]
	x1 <- xmat[, match(cissnp, bim$V2)]
	x <- xmat[, match(tsnps, bim$V2)]
	y <- resphen[, match(pn, colnames(resphen))]

	temp1 <- subset(s, select=c(snp1, chr1))
	temp2 <- subset(s, select=c(snp2, chr2))
	names(temp1) <- names(temp2) <- c("snp", "chr")
	chrkey <- rbind(temp1, temp2)
	chrkey <- subset(chrkey, !duplicated(snp))

	l <- list()
	for(i in 1:ncol(x))
	{
		print(table(x1, x[,i]))
		gp <- tapply(y, list(x1, x[,i]), function(x) mean(x, na.rm=T))
		gp <- gp - min(gp, na.rm=T)
		print(gp)
		title <- paste(s$probegene[1], "chr", subset(chrkey, snp == cissnp)$chr, "x", subset(chrkey, snp == tsnps[i])$chr)
		l[[i]] <- plot3dGp(gp, title, cissnp, tsnps[i], z)
	}
	do.call(grid.arrange, l)
}


load("~/repo/eQTL-2D/analysis/replication2_summary.RData")


