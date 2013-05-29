library(lattice)
library(latticeExtra)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(plyr)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL)
{
	require(grid)

	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)

	numPlots = length(plots)

	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) 
	{
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
		ncol = cols, nrow = ceiling(numPlots/cols))
	}

	if (numPlots==1) 
	{
		print(plots[[1]])
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

		# Make each plot, in the correct location
		for (i in 1:numPlots)
		{
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
			layout.pos.col = matchidx$col))
		}
	}
}


plot3d <- function(res, geno, phen, z=-45)
{
	a <- tapply(phen[,res$probeid], list(geno[,res$pos1], geno[, res$pos2]), function(x) mean(x, na.rm=T))
	a <- a - min(a, na.rm=T)
	b <- table(geno[,res$pos1], geno[, res$pos2])
	print(b)
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

plot3dProbe <- function(res, pg, geno, phen, z=45)
{
	a <- subset(res, probegene == pg)
	l <- list()
	for(i in 1:nrow(a))
	{
		l[[i]] <- plot3d(a[i,], xmat, resphen, z)[[3]]
	}
	do.call(grid.arrange, l)
}

subSigPg <- function(pg, s=sig)
{
	return(subset(s, probegene==pg))
}

# A better way to plot these 3D graphs

plot3dHairballs <- function(sig, pg, cissnp, resphen, xmat, bim, z=45)
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


#=================================================================================================#
#=================================================================================================#


load("~/repo/eQTL-2D/analysis/interaction_list_replication_summary.RData")
load("~/repo/eQTL-2D/data/residuals_all.RData")
load("~/repo/eQTL-2D/data/clean_geno_final.RData")


pdf("~/repo/eQTL-2D/analysis/images/TMEM149.pdf", width=20, height=20)
plot3dHairballs(sig, "TMEM149", "rs8106959", resphen, xmat, bim, z=-45)
dev.off()

pdf("~/repo/eQTL-2D/analysis/images/MBNL1.pdf", width=20, height=20)
plot3dHairballs(sig, "MBNL1", "rs13069559", resphen, xmat, bim, z=-45)
dev.off()

pdf("~/repo/eQTL-2D/analysis/images/CAST.pdf", width=20, height=20)
plot3dHairballs(sig, "CAST", "rs7733671", resphen, xmat, bim, z=45)
dev.off()

pdf("~/repo/eQTL-2D/analysis/images/TRAPPC5.pdf", width=20, height=20)
plot3dHairballs(sig, "TRAPPC5", "rs17159840", resphen, xmat, bim, z=45)
dev.off()

pdf("~/repo/eQTL-2D/analysis/images/HMBOX1.pdf", width=15, height=15)
plot3dProbe(sig, "HMBOX1", xmat, resphen)
dev.off()

pdf("~/repo/eQTL-2D/analysis/images/NAPRT1.pdf", width=15, height=15)
plot3dHairballs(sig, "NAPRT1", "rs2123758", resphen, xmat, bim, z=-45)
dev.off()

pdf("~/repo/eQTL-2D/analysis/images/ADK.pdf", width=15, height=15)
plot3dProbe(sig, "ADK", xmat, resphen)
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



#=================================================================================================#
#=================================================================================================#

sigb <- subset(sig, pnest_egcut > -log10(0.05/380) | pnest_fehr > -log10(0.05/380))

# Let's look at ADK

pdf("~/repo/eQTL-2D/analysis/images/ADK.pdf", width=15, height=15)
plot3dProbe(sig, "ADK", xmat, resphen, z=45)
dev.off()

# both SNPs close together on chromosome 10
# first in an intron of a 5' gene of ADK called VCL
# Second in the last intron of ADK

# ADK = adenosine kinase, expressed everywhere
# VCL = Vinculin - expressed in muscle tissue

# VCL snp is probably modifying promotor activity for ADK


load("~/repo/eQTL-2D/data/probeinfo_all.RData")
load("~/repo/eQTL-2D/filtering/marginal_lists/marginal_list.RData")
temp <- subset(probeinfo_all, select=c(PROBE_ID, ILMN_GENE))

marginal_list <- merge(marginal_list, temp, by.x="probename", by.y="PROBE_ID", all.x=T)
subset(marginal_list, snp %in% c("rs2395095", "rs10824092"))

# Both have an additive effect on ADK




#=================================================================================================#
#=================================================================================================#

# CTSC
plot3dProbe(sig, "CTSC", xmat, resphen, z=135)






#=================================================================================================#
#=================================================================================================#



# Distribution of additive vs non-additive variance

a <- as.data.frame(do.call(rbind, sig$vc))
sig$varA <- a$.a + a$a.
sig$varD <- a$.d + a$d.
sig$varI <- a$aa + a$ad + a$da + a$dd
prop <- subset(sig, select=c(varA, varI, varD))
prop <- prop[order(prop$varA+prop$varI+prop$varD, decreasing=T),]
prop$index <- 1:nrow(prop)
prop <- melt(prop, id=c("index"))
prop <- prop[nrow(prop):1, ]
prop$variable <- factor(prop$variable, levels=c("varI", "varD", "varA"))
levels(prop$variable) <- c("Interaction", "Dominance", "Additive")
p1 <- ggplot(prop, aes(y=value, x=index)) + 
	geom_bar(stat="identity", aes(fill=variable, colour=variable), position=position_stack(width=0)) +
	scale_fill_brewer("Variance component") +
	scale_colour_brewer("Variance component") +
	ylab("Phenotypic variance") + xlab("") +
	coord_flip() + theme(legend.position = "none") +
	theme(axis.text.y=element_text(size=0), axis.ticks.y=element_line(size=0))
ggsave("~/repo/eQTL-2D/analysis/images/proportion_additive.pdf", width=10, height=10)


a <- as.data.frame(do.call(rbind, sig$vc))
sig$varA <- a$.a + a$a.
sig$varD <- a$.d + a$d.
sig$varI <- a$aa + a$ad + a$da + a$dd
sig$varG <- with(sig, varA + varD + varI)
sig$varA <- sig$varA / sig$varG
sig$varD <- sig$varD / sig$varG
sig$varI <- sig$varI / sig$varG
prop <- subset(sig, select=c(varA, varI, varD))
prop <- prop[order(prop$varA, decreasing=T),]
prop$index <- 1:nrow(prop)
prop <- melt(prop, id=c("index"))
prop <- prop[nrow(prop):1, ]
prop$variable <- factor(prop$variable, levels=c("varI", "varD", "varA"))
levels(prop$variable) <- c("Interaction", "Dominance", "Additive")
p2 <- ggplot(prop, aes(y=value, x=index)) + 
	geom_bar(stat="identity", aes(fill=variable, colour=variable), position=position_stack(width=0)) +
	scale_fill_brewer("Variance component") +
	scale_colour_brewer("Variance component") +
	ylab("Proportion of genetic variance") + xlab("") +
	coord_flip() + 
	theme(legend.justification=c(1,0), legend.position=c(1,0)) + 
	theme(axis.text.y=element_text(size=0), axis.ticks.y=element_line(size=0))
ggsave("~/repo/eQTL-2D/analysis/images/proportion_genetic.pdf", width=10, height=10)

pdf(file="~/repo/eQTL-2D/analysis/images/variance_components.pdf", width=10, height=10)
multiplot(p1, p2, cols=2)
dev.off()

