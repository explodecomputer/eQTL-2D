library(lattice)
library(latticeExtra)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(plyr)
library(xtable)


#=================================================================================================#
#=================================================================================================#


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
		g$s1 <- a$snp1[i]
		g$s2 <- a$snp2[i]
		g$chr1 <- a$chr1[i]
		g$chr2 <- a$chr2[i]
		g$index <- i
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


plot3dGp <- function(gp, title="", snp1="SNP1", snp2="SNP2", z=-45)
{
	p <- cloud(
		gp, # This has to be a table!?
		panel.3d.cloud=panel.3dbars, 
		col="black", 
		col.facet=c("#e5f5e0", "#A1D99B", "#31A354"), 
		# col.facet=c("#31A354"), 
		xbase=0.4, 
		ybase=0.4,
		xlab=list(
			label=snp1, 
			cex=0.8), 
		ylab=list(
			label=snp2, 
			cex=0.8), 
		scales=list(
			arrows=F, 
			z = list(cex = 0), 
			y = list(cex = 0.5), 
			x=list(cex=0.5)),
		zlab="",
		cex.title = 0.5,
		screen = list(
			z = z, 
			x = -60, 
			y = 3),
		main = list(
			label=title, 
			cex=0.8)
	)
	return(p)
}


plot3dGpGrid <- function(dat, dataset, cischr, z)
{
	l <- list()
	index <- unique(dat$index)
	j <- 1
	for(i in index)
	{
		a <- subset(dat, index == i & cohort == dataset)
		print(head(a))
		gp <- matrix(a$y, 3, 3)
		gp <- gp - min(gp, na.rm=T)
		gp <- as.table(gp)
		rownames(gp) <- colnames(gp) <- c("0", "1", "2")
		m <- paste("Chr", a$chr1[1], "x", a$chr2[1], "\n", a$s1[1], "x", a$s2[1])
		sl1 <- ifelse(a$chr1[1] == cischr, "cis", "trans")
		sl2 <- ifelse(a$chr2[1] == cischr, "cis", "trans")
		zi <- ifelse(is.na(z[i]), -45, z[i])
		l[[i]] <- plot3dGp(as.table(gp), title=m, snp1=sl1, snp2=sl2, z=zi)
	}
	do.call(grid.arrange, l)
}


#=================================================================================================#
#=================================================================================================#


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


#=================================================================================================#
#=================================================================================================#


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


#=================================================================================================#
#=================================================================================================#


MBNL1 <- subset(datHeatmapGp(subset(sig, probegene == "MBNL1")), cohort == "BSGS")
pdf(file="~/repo/eQTL-2D/analysis/images/MBNL1_3d.pdf", width=10, height=10)
plot3dGpGrid(MBNL1, "BSGS", 3, c(45))
dev.off()

TMEM149 <- subset(datHeatmapGp(subset(sig, probegene == "TMEM149")), cohort == "BSGS")
pdf(file="~/repo/eQTL-2D/analysis/images/TMEM149_3d.pdf", width=10, height=12.5)
plot3dGpGrid(TMEM149, "BSGS", 19, c(-45, rep(135, 14), 45, 135, 135))
dev.off()

CAST <- subset(datHeatmapGp(subset(sig, probegene == "CAST")), cohort == "BSGS")
pdf(file="~/repo/eQTL-2D/analysis/images/CAST_3d.pdf", width=10, height=10)
plot3dGpGrid(CAST, "BSGS", 5, c(rep(135, 4), 225, 135, 135, 225, 225, 135, 135, 135, 225, 135, 135))
dev.off()

NAPRT1 <- subset(datHeatmapGp(subset(sig, probegene == "NAPRT1")), cohort == "BSGS")
pdf(file="~/repo/eQTL-2D/analysis/images/NAPRT1_3d.pdf", width=7.5, height=7.5)
plot3dGpGrid(NAPRT1, "BSGS", 8, c(135, 135, 135, 135, 45, 135, 135, 45))
dev.off()

TRAPPC5 <- subset(datHeatmapGp(subset(sig, probegene == "TRAPPC5")), cohort == "BSGS")
pdf(file="~/repo/eQTL-2D/analysis/images/TRAPPC5_3d.pdf", width=10, height=10)
plot3dGpGrid(TRAPPC5, "BSGS", 19, c(-45, 225, -45, -45, -45, 225))
dev.off()


#=================================================================================================#
#=================================================================================================#

