library(lattice)
library(latticeExtra)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(plyr)


plot3d <- function(res, geno, phen, z=45)
{
	a <- tapply(phen[,res$probeid], list(geno[,res$pos1], geno[, res$pos2]), mean)
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











