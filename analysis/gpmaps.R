library(lattice)
library(latticeExtra)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(plyr)
library(xtable)


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
		g <- GpMod(g, i, "BSGS", rotateGp, rotateGp)
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
	p <- ggplot(l, aes(x=snp2, y=snp1, fill=y)) + 
	geom_tile() + 
	facet_grid(code3 ~ cohort) + 
	theme(strip.text = element_text(size = 6),
		axis.text = element_text(size = 6),
		axis.title = element_text(size = 8), 
		legend.position="none") + 
	xlab("SNP 2") + 
	ylab("SNP 1")
	return(p)
}

plotTileGp2 <- function(l)
{
	p <- ggplot(l, aes(x=snp2, y=snp1, fill=y)) + 
	geom_tile() + 
	facet_grid(code3 ~ cohort) + 
	theme(strip.text = element_text(size = 6), 
		axis.ticks.y = element_line(0),
		axis.text.y = element_text(size=0),
		legend.position = "none") + 
	xlab("SNP 2") + 
	ylab(NULL)
	return(p)
}


convertToScientific <- function(x)
{

}


load("~/repo/eQTL-2D/analysis/interaction_list_meta_analysis.RData")
sig <- subset(meta, filter !=3)

# Results table

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



# Plots

l <- datHeatmapGp(bsig)
l <- GpMod(l, 2, "EGCUT", rotateGp, rotateGp, rotateGp)
l <- GpMod(l, 3, "EGCUT", t)
l <- GpMod(l, 5, "BSGS", flipGpCol)
l <- GpMod(l, 6, "EGCUT", flipGpRow)
l <- GpMod(l, 10, "Fehrmann", flipGpCol)
l <- GpMod(l, 11, "EGCUT", flipGpRow)
l <- GpMod(l, 12, "EGCUT", flipGpRow) 
l <- GpMod(l, 13, "EGCUT", flipGpRow)
l <- GpMod(l, 15, "EGCUT", flipGpRow)
l <- GpMod(l, 22, "BSGS", flipGpRow)
l <- GpMod(l, 22, "EGCUT", flipGpRow)
l <- GpMod(l, 22, "Fehrmann", flipGpRow)


pdf(file="~/repo/eQTL-2D/analysis/images/gpBonfRep.pdf", width=4, height=7)
plot1 <- plotTileGp(subset(l, code3 %in% 1:15))
plot2 <- plotTileGp(subset(l, code3 %in% 16:30))
grid.arrange(plot1, plot2, ncol=2)
dev.off()
plot1


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



m <- subset(marginal_list, pval < 10^-15.51 & chr %in% 1:22)
dim(m)
m1 <- subset(m, !duplicated(probename) & probename %in% probeinfo$PROBE_ID)
dim(m1)
m2 <- merge(m1, probeinfo, by.x="probename", by.y="PROBE_ID")
m3 <- subset(m2, CHR %in% 1:22)
dim(m3)
m3 <- subset(m3, !duplicated(ILMN_GENE))
dim(m3)

517 to 453

vg <- 







