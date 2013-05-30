library(lattice)
library(latticeExtra)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(plyr)


plot3dGp <- function(gp, title="", snp1="SNP1", snp2="SNP2", z=-45)
{
	p <- cloud(
		gp, 
		panel.3d.cloud=panel.3dbars, 
		col="black", 
		col.facet=c("#e5f5e0", "#A1D99B", "#31A354"), 
		xbase=0.6, ybase=0.6,
		xlab=list(label=snp1, cex=0.8), ylab=list(label=snp2, cex=0.8), zlab="",
		default.scales=list(arrows=F, z = list(cex = 0), y = list(cex = 0.5), x=list(cex=0.5)),
		cex.title = 0.5,
		screen = list(z = z, x = -60, y = 3),
		main = title
	)
}


plot3dHairballsRep <- function(sig, pg, cissnp, cohort, z=45)
{
	# make table of mean phenotypes
	# AA cis vs 1:nalleles trans
	# Aa cis vs 1:nalleles trans
	# aa cis vs 1:nalleles trans

	s <- subset(sig, probegene==pg)
	snps <- unique(c(s$snp1, s$snp2))
	stopifnot(cissnp %in% snps)
	tsnps <- snps[! snps %in% cissnp]

	# transpose the GP map for all where the cis snp is snp2

	temp1 <- subset(s, select=c(snp1, chr1))
	temp2 <- subset(s, select=c(snp2, chr2))
	names(temp1) <- names(temp2) <- c("snp", "chr")
	chrkey <- rbind(temp1, temp2)
	chrkey <- subset(chrkey, !duplicated(snp))

	nom <- paste("gcm", cohort, sep="")

	l <- list()
	for(i in 1:nrow(s))
	{
		print(gp <- s[[nom]][[i]])
		snp1 <- cissnp
		snp2 <- s$snp2[i]

		if(s$snp1[i] != cissnp)
		{
			gp <- t(gp)
			snp2 <- s$snp1[i]
		}
		gp <- gp - min(gp, na.rm=T)

		title <- paste("chr", subset(chrkey, snp == cissnp)$chr, "x", subset(chrkey, snp == snp2)$chr)
		l[[i]] <- plot3dGp(gp, title, cissnp, snp2, z)
	}
	do.call(grid.arrange, l)
}

normaliseGp <- function(g)
{
	a <- ddply(g, .(cohort), function(x)
	{
		x <- mutate(x)
		x$y <- x$y - min(x$y, na.rm=T)
		x$y <- x$y / max(x$y, na.rm=T)			
		return(x)
	})
	return(a)
}

datHeatmapGp <- function(sig)
{
	a <- subset(sig, select=c(snp1, snp2, chr1, chr2, probegene, probename, gcm, gcm_fehr, gcm_egcut, code))

	l <- list()
	for(i in 1:nrow(a))
	{
		g1 <- expand.grid(snp1=0:2, snp2=0:2)
		g1$y <- c(a$gcm[[i]])
		g1$cohort <- "BSGS"
		g2 <- expand.grid(snp1=0:2, snp2=0:2)
		g2$y <- c(a$gcm_egcut[[i]])
		g2$cohort <- "EGCUT"
		g3 <- expand.grid(snp1=0:2, snp2=0:2)
		g3$y <- c(a$gcm_fehr[[i]])
		g3$cohort <- "Fehrmann"
		g <- normaliseGp(rbind(g1, g2, g3))

		g$code <- a$code[i]
		g$code2 <- a$probegene[i]
		g$code3 <- i
		l[[i]] <- g
	}
	l <- rbind.fill(l)
	return(l)
}

GpMod <- function(l, code, cohort, ...)
{
	index <- l$code3 == code & l$cohort == cohort
	mat <- matrix(l$y[index], 3, 3, byrow=FALSE)
	print(mat)
	print(l$y[index])
	funcs <- list(...)
	for(i in 1:length(funcs))
	{
		mat <- funcs[[i]](mat)
	}
	print(mat)
	l$y[index] <- c((mat))
	print(l$y[index])
	return(l)
}

flipGpRow <- function(mat)
{
	return(mat[3:1, ])
}

flipGpCol <- function(mat)
{
	return(mat[, 3:1])
}

rotateGp <- function(mat)
{
	return(t(mat[3:1, ]))
}

plotTileGp <- function(l)
{
	p <- ggplot(l, aes(x=snp2, y=snp1, fill=y)) + geom_tile() + facet_grid(cohort ~ code3) + theme(strip.text = element_text(size = 6))
	return(p)
}


load("~/repo/eQTL-2D/analysis/interaction_list_replication_summary.RData")

l <- datHeatmapGp(subset(sig, pnest_fehr > -log10(0.05/500) & pnest_egcut > -log10(0.05/500)))
l <- GpMod(l, 1, "BSGS", flipGpCol)
l <- GpMod(l, 1, "Fehrmann", flipGpRow)
l <- GpMod(l, 10, "BSGS", flipGpRow, rotateGp, rotateGp, rotateGp)
l <- GpMod(l, 11, "BSGS", flipGpRow)
l <- GpMod(l, 12, "BSGS", rotateGp, rotateGp)
l <- GpMod(l, 13, "BSGS", rotateGp, rotateGp)
l <- GpMod(l, 13, "Fehrmann", flipGpCol)
l <- GpMod(l, 14, "BSGS", rotateGp, rotateGp)
l <- GpMod(l, 15, "BSGS", rotateGp, rotateGp)
l <- GpMod(l, 16, "BSGS", rotateGp, rotateGp)
l <- GpMod(l, 17, "BSGS", rotateGp, rotateGp)
l <- GpMod(l, 2, "Fehrmann", rotateGp)
l <- GpMod(l, 2, "BSGS", rotateGp, rotateGp, rotateGp)
l <- GpMod(l, 3, "BSGS", rotateGp, rotateGp)
l <- GpMod(l, 4, "BSGS", rotateGp, rotateGp)
l <- GpMod(l, 5, "BSGS", rotateGp, rotateGp)
l <- GpMod(l, 6, "BSGS", t)
l <- GpMod(l, 7, "BSGS", rotateGp, rotateGp)
l <- GpMod(l, 9, "BSGS", rotateGp, rotateGp)

pdf(file="~/repo/eQTL-2D/analysis/images/gpBonfRep.pdf", width=20, height=4)
plotTileGp(l)
dev.off()

temp <- subset(sig, probegene == "MBNL1")
temp <- subset(temp, !duplicated(paste(chr1, chr2)) & !is.na(pnest_egcut) & !is.na(pnest_fehr))
mbnl1 <- datHeatmapGp(temp)
mbnl1 <- GpMod(mbnl1, 1, "BSGS", flipGpCol)
mbnl1 <- GpMod(mbnl1, 2, "BSGS", rotateGp, rotateGp)
mbnl1 <- GpMod(mbnl1, 3, "BSGS", flipGpCol)
mbnl1 <- GpMod(mbnl1, 4, "BSGS", rotateGp, rotateGp)
mbnl1 <- GpMod(mbnl1, 5, "BSGS", flipGpCol)
mbnl1 <- GpMod(mbnl1, 6, "BSGS", flipGpRow)
mbnl1 <- GpMod(mbnl1, 6, "Fehrmann", flipGpCol)
mbnl1 <- GpMod(mbnl1, 7, "BSGS", rotateGp, rotateGp)
mbnl1 <- GpMod(mbnl1, 8, "BSGS", rotateGp, rotateGp)
mbnl1 <- GpMod(mbnl1, 9, "BSGS", flipGpCol)
mbnl1 <- GpMod(mbnl1, 10, "BSGS", flipGpCol)
mbnl1 <- GpMod(mbnl1, 11, "BSGS", flipGpCol)
mbnl1 <- GpMod(mbnl1, 12, "BSGS", rotateGp, rotateGp)
mbnl1 <- GpMod(mbnl1, 12, "Fehrmann", flipGpRow)
mbnl1 <- GpMod(mbnl1, 13, "BSGS", flipGpCol)
mbnl1 <- GpMod(mbnl1, 14, "BSGS", flipGpCol)
mbnl1 <- GpMod(mbnl1, 14, "Fehrmann", flipGpRow)

index <- match(mbnl1$code, sig$code)
replicates <- rep("", length(index))
replicates[with(sig[index, ], pnest_fehr > upper_fehr | pnest_egcut > upper_egcut)] <- "*"
replicates[with(sig[index, ], pnest_fehr > upper_fehr & pnest_egcut > upper_egcut)] <- "**"
mbnl1$code4 <- paste(sig$chr1[index], "x", sig$chr2[index], replicates)
mbnl1 <- subset(mbnl1, !is.na(y))
mbnl1$code3 <- mbnl1$code4

plotTileGp(mbnl1)



plot3dHairballsRep(sig, "TMEM149", "rs8106959", "", z=-45)
plot3dHairballsRep(sig, "TMEM149", "rs8106959", "_egcut", z=45)
plot3dHairballsRep(sig, "TMEM149", "rs8106959", "_fehr", z=45)

pdf("~/repo/eQTL-2D/analysis/images/MBNL1.pdf", width=20, height=20)
plot3dHairballsRep(subset(sig, code %in% mbnl1$code), "MBNL1", "rs13069559", "", z=-135)
dev.off()

plot3dHairballsRep(subset(sig, code %in% mbnl1$code), "MBNL1", "rs13069559", "_egcut", z=45)
plot3dHairballsRep(subset(sig, code %in% mbnl1$code), "MBNL1", "rs13069559", "_fehr", z=45)

adk <- subset(sig, probegene=="ADK")[1,]

par(mfrow=c(2,2))
print(plot3dGp(adk$gcm[[1]], z=45))
print(plot3dGp(adk$gcm_fehr[[1]], z=45))
print(plot3dGp(adk$gcm_egcut[[1]], z=45))

