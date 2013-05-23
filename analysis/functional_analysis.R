library(lattice)
library(latticeExtra)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(plyr)


plot3d <- function(res, geno, phen, z=45)
{
	a <- tapply(phen[,res$probeid], list(geno[,res$pos1], geno[, res$pos2]), function(x) mean(x, na.rm=T))
	b <- table(geno[,res$pos1], geno[, res$pos2])
	p <- cloud(
		a, 
		panel.3d.cloud=panel.3dbars, 
		col="black", 
		col.facet=c("#e5f5e0", "#A1D99B", "#31A354"), 
		xbase=0.6, ybase=0.6,
		xlab="SNP1", ylab="SNP2", zlab="y",
		default.scales=list(arrows=F),
		screen = list(z = z, x = -60, y = 3),
		main = paste(res$chr1, res$chr2)
	)
	return(list(a,b,p))
}


plot3dProbe <- function(res, pg, geno, phen, z=45)
{
	a <- subset(res, probegene == pg)
	l <- list()
	for(i in 1:nrow(a))
	{
		l[[i]] <- plot3d(a[i,], xmat, resphen)[[3]]
	}
	do.call(grid.arrange, l)
}


#=================================================================================================#
#=================================================================================================#


load("~/repo/eQTL-2D/analysis/replication_summary.RData")
load("~/repo/eQTL-2D/data/residuals_all.RData")
load("~/repo/eQTL-2D/data/clean_geno_final.RData")


pdf("~/repo/eQTL-2D/analysis/images/TMEM149.pdf", width=15, height=15)
plot3dProbe(sig, "TMEM149", xmat, resphen)
dev.off()

pdf("~/repo/eQTL-2D/analysis/images/MBNL1.pdf", width=15, height=15)
plot3dProbe(sig, "MBNL1", xmat, resphen)
dev.off()

pdf("~/repo/eQTL-2D/analysis/images/CAST.pdf", width=15, height=15)
plot3dProbe(sig, "CAST", xmat, resphen)
dev.off()

pdf("~/repo/eQTL-2D/analysis/images/TRAPPC5.pdf", width=15, height=15)
plot3dProbe(sig, "TRAPPC5", xmat, resphen)
dev.off()

pdf("~/repo/eQTL-2D/analysis/images/HMBOX1.pdf", width=15, height=15)
plot3dProbe(sig, "HMBOX1", xmat, resphen)
dev.off()

pdf("~/repo/eQTL-2D/analysis/images/NAPRT1.pdf", width=15, height=15)
plot3dProbe(sig, "NAPRT1", xmat, resphen)
dev.off()


#=================================================================================================#
#=================================================================================================#

# Get list of genes for each SNP
# Use HaploReg to get the gene list
# Use Partner Hunter or Set Distiller in GeneDecks to get the list of associated genes

# Make SNP list for HaploReg

makeHrSnpList <- function(sig, filename)
{
	hr <- unique(c(sig$snp1, sig$snp2))
	write.table(hr, file=filename, row=F, col=F, qu=F)
}

readHrResults <- function(filename)
{
	hr2 <- scan(filename, what=character(), na.strings=".")
	h <- hr2[1:31]
	hr2 <- hr2[-c(1:31)]
	hr2 <- as.data.frame(matrix(hr2, ncol=31, byrow=T), stringsAsFactors=FALSE)
	names(hr2) <- h
	hr2$r2 <- as.numeric(hr2$r2)
	hr2$chr <- as.numeric(hr2$chr)
	hr2$pos <- as.numeric(hr2$pos)
	hr2$is_query_snp <- as.numeric(hr2$is_query_snp)*-1
	hr2 <- hr2[order(hr2$chr, hr2$pos), ]
	return(hr2)
}

getRsGene <- function(hr, sig)
{
	a <- subset(hr2, is_query_snp == 1)
	a$GENCODE_name
	a <- subset(a, select=c(rsID, GENCODE_name))
	head(sig)
	sig <- merge(sig, a, by.x="snp1", by.y="rsID", all.x=T)
	sig <- merge(sig, a, by.x="snp2", by.y="rsID", all.x=T, suff=c(1,2))
	return(sig)	
}

getProbeGeneList <- function(sig, pg)
{
	a <- subset(sig, probegene==pg)
	print(with(a, table(c(probegene, GENCODE_name1, GENCODE_name2))))
	genelist <- with(a, unique(c(probegene, GENCODE_name1, GENCODE_name2)))
	for(i in 1:length(genelist)) 
	{
		cat(genelist[i], "\n")
	}
}

hr2 <- readHrResults("HaploReg_results.txt")
mbnl1 <- subset(sig, probegene == "MBNL1")
ms <- unique(c(mbnl1$snp1, mbnl1$snp2))
a <- subset(hr2, rsID %in% ms)
plot(log(diff(hr2$pos)[diff(hr2$pos)>0]))
getProbeGeneList(sig, "TMEM149")
getProbeGeneList(sig, "MBNL1")
getProbeGeneList(sig, "CAST")
getProbeGeneList(sig, "TRAPPC5")
getProbeGeneList(sig, "NAPRT1")


# Nothing very significant shows up. Kostya will look at enrichment of motifs



#=================================================================================================#
#=================================================================================================#


# Does the number of trans alleles have an effect?

cumTransAlleles <- function(sig, pg, cissnp, resphen, xmat, bim)
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
	x_code <- x
	info <- data.frame(tsnps, ref=2)

	for(i in 1:ncol(x))
	{
		print(table(x1, x[,i]))
		gp <- tapply(y, list(x1, x[,i]), function(x) mean(x, na.rm=T))
		print(gp)
		print(info$ref[i] <- ifelse(gp[1,3] > gp[3,3], 0, 2))
		# index <- x_code[,i] == info$ref[i]
		# x_code[index,i] <- 1
		# x_code[!index,i] <- 0
		if(info$ref[i] == 0)
		{
			x_code[,i] <- -1*(x_code[,1]-1)+1
		}

	}

	x_code <- apply(x_code, 1, sum)
	l <- list()
	l$gp <- tapply(y, list(x1, x_code), function(x) mean(x, na.rm=T))
	l$count <- table(list(x1, x_code))
	l$se <- tapply(y, list(x1, x_code), function(x) sd(x, na.rm=T)) / sqrt(l$count)
	return(l)
}


cta <- cumTransAlleles(sig, "TMEM149", "rs8106959", resphen, xmat, bim)
a <- melt(cta$gp)
b <- melt(cta$se)
a$se <- b$value
ggplot(a, aes(x=Var2, y=value)) + geom_line() + facet_grid(Var1 ~ .) + geom_errorbar(aes(ymax=value+se, ymin=value-se))














