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
		xlab=snp1, ylab=snp2, zlab="y",
		default.scales=list(arrows=F),
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

		title <- paste(s$probegene[1], "chr", subset(chrkey, snp == cissnp)$chr, "x", subset(chrkey, snp == snp2)$chr)
		l[[i]] <- plot3dGp(gp, title, cissnp, snp2, z)
	}
	do.call(grid.arrange, l)
}

rotateGp <- function(mat)
{
	if(!is.na(mat[1]))
	{
		return(mat[3:1, 3:1])
	} else {
		return(NA)
	}
}

normaliseGp <- function(g)
{
	a <- ddply(g, .(cohort), function(x)
	{
		x <- mutate(x)
		if(length(x) == 1 & is.na(x[1]))
		{
			return(x)
		} else {
			x$y <- x$y - min(x$y, na.rm=T)
			x$y <- x$y / max(x$y, na.rm=T)			
		}
		return(x)
	})
	print(head(a))
	return(a)
}

plotHeatmapGp <- function(sig)
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
		g$code2 <- paste(i, a$probegene[i])
		l[[i]] <- g
	}
	l <- rbind.fill(l)
	print(head(l))
	p <- ggplot(l, aes(x=snp1, y=snp2, fill=y)) + geom_tile() + facet_grid(cohort ~ code2) + theme(strip.text = element_text(size = 6))

	return(p)
}


load("~/repo/eQTL-2D/analysis/interaction_list_replication_summary.RData")

plot3dHairballsRep(sig, "TMEM149", "rs8106959", "", z=45)
plot3dHairballsRep(sig, "TMEM149", "rs8106959", "_egcut", z=45)
plot3dHairballsRep(sig, "TMEM149", "rs8106959", "_fehr", z=45)

plot3dHairballsRep(sig, "MBNL1", "rs13069559", "", z=-45)
plot3dHairballsRep(sig, "MBNL1", "rs13069559", "_egcut", z=-45)
plot3dHairballsRep(sig, "MBNL1", "rs13069559", "_fehr", z=-45)

adk <- subset(sig, probegene=="ADK")[1,]

par(mfrow=c(2,2))
print(plot3dGp(adk$gcm[[1]], z=45))
print(plot3dGp(adk$gcm_fehr[[1]], z=45))
print(plot3dGp(adk$gcm_egcut[[1]], z=45))


pdf(file="~/repo/eQTL-2D/analysis/images/gpBonfRep.pdf", width=20, height=4)
plotHeatmapGp(subset(sig, pnest_fehr > -log10(0.05/500) & pnest_egcut > -log10(0.05/500)))
dev.off()
