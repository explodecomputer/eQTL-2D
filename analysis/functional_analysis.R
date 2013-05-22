library(lattice)
library(latticeExtra)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(plyr)


plot3D <- function(res, geno, phen, z=45)
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






load("~/repo/eQTL-2D/analysis/replication_summary.RData")

a <- subset(sig, probegene == "TMEM149")

l <- list()
for(i in 1:nrow(a))
{
	l[[i]] <- plot3D(a[i,], xmat, resphen)[[3]]
}


b <- subset(sig, probegene == "MBNL1")

l2 <- list()
for(i in 1:nrow(b))
{
	l2[[i]] <- plot3D(a[i,], xmat, resphen)[[3]]
}

pdf("~/repo/eQTL-2D/analysis/TMEM149.pdf", width=15, height=15)
do.call(grid.arrange, l)
dev.off()

pdf("~/repo/eQTL-2D/analysis/MBNL1.pdf", width=15, height=15)
do.call(grid.arrange, l2)
dev.off()




