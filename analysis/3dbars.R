library(lattice)
library(latticeExtra)
library(noia)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(plyr)

# Load in data

load("~/repo/eQTL-2D/data/residuals_all.RData")
load("~/repo/eQTL-2D/filtering/filtered_by_chr/sig_all_probes.RData")
load("~/repo/eQTL-2D/data/clean_geno_final.RData")



#' What is the biggest variance component for each interaction?
#'
#' For each interaction decompose the variance using NOIA and 
#' find the epistatic VC with the largest variance
#' @param sig A single row from sig
#' @param geno Genotype matrix
#' @param phen Residuals for each probe
#' @param loud=FALSE Print full linear model?
#' @export
#' @return aa, ad, da or dd
#' @alias
#' @examples \dontrun{
#' varianceComponents(sig[1,], xmat, resphen, TRUE)
#'}
varianceComponents <- function(sig, geno, phen, loud=FALSE)
{
	a <- linearRegression(
		phen = phen[,sig$probeid],
		gen  = cbind(geno[,sig$pos1], geno[, sig$pos2]) + 1
	)
	if(loud) print(a)
	vars <- a$variances
	vars <- vars[names(vars) %in% c("aa", "ad", "da", "dd")]

	nom <- names(vars)[which.max(vars)]

	nom <- unlist(strsplit(toupper(nom), split=""))
	nom <- paste(nom, collapse=" x ")

	return(nom)
}

#' What is the biggest variance component for each interaction?
#'
#' For each interaction decompose the variance using NOIA and 
#' find the epistatic VC with the largest variance
#' @param sig A single row from sig
#' @param geno Genotype matrix
#' @param phen Residuals for each probe
#' @param loud=FALSE Print full linear model?
#' @export
#' @return aa, ad, da or dd
#' @alias
#' @examples \dontrun{
#' varianceComponentsBreakdown(sig[1,], xmat, resphen, TRUE)
#'}
varianceComponentsBreakdown <- function(sig, geno, phen, loud=FALSE)
{
	a <- linearRegression(
		phen = phen[,sig$probeid],
		gen  = cbind(geno[,sig$pos1], geno[, sig$pos2]) + 1
	)
	if(loud) print(a)
	vars <- as.data.frame(as.list(a$variances))[,-1]
	sig <- cbind(sig, vars)

	return(sig)
}

#' Run \link{varianceComponents} on all interaction pairs
#'
#' @param sig A single row from sig
#' @param geno Genotype matrix
#' @param phen Residuals for each probe
#' @export
#' @return sig with an extra column called vc
#' @examples \dontrun{
#' sig <- getVc(sig, xmat, resphen)
#' a <- table(sig$vc)
#' mat <- rbind(a, rep(mean(a),4))
#' chisq.test(mat, simulate.p.value = TRUE, B = 10000)
#'}
getVc <- function(sig, geno, phen)
{
	sig$vc <- NA
	for(i in 1:nrow(sig))
	{
		cat(i, "of", nrow(sig), "\n")
		sig$vc[i] <- varianceComponents(sig[i,], xmat, resphen)
	}
	return(sig)
}

getVcBreakdown <- function(sig, geno, phen)
{
	l <- list()
	for(i in 1:nrow(sig))
	{
		cat(i, "of", nrow(sig), "\n")
		l[[i]] <- varianceComponentsBreakdown(sig[i,], xmat, resphen)
	}
	l <- rbind.fill(l)
	return(l)
}

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


qqPlot <- function(pval)
{
	o <- -log10(sort(pval,decreasing=F))
	e <- -log10(ppoints(length(pval)))
	ret <- data.frame(expected=e, observed=o)
	print(plot(NULL, type="n", ylim=c(0, max(c(o,e))), xlim=c(0, max(c(o,e)))))
	print(plot(ret))
	print(abline(a=0, b=1))
	return(ret)
}

a <- with(newsig, qqPlot(10^-replication_pnest))
plot(observed ~ expected, data=a[-c(1:20),])
abline(a=0, b=1)
maf <- apply(xmat, 2, function(x) {
	sum(x, na.rm=T) / (2*length(x))
})
sig$maf1 <- maf[sig$pos1]
sig$maf2 <- maf[sig$pos2]
sig$mmaf <- (sig$maf1 + sig$maf2) / 2
sig$propG[sig$propG >= 1] <- 0.95
sig$varA <- sig$propG * sig$propA
sig$varI <- sig$propG * (1-sig$propA)
sig <- getVc(sig, xmat, resphen)
sig <- getVcBreakdown(sig, xmat, resphen)
save(sig, file="~/repo/eQTL-2D/analysis/sig_analysis.RData")
# load("~/repo/eQTL-2D/analysis/sig_analysis.RData")


############
# 3D plots #
############

sig2 <- sig[order(sig$propA, decreasing=FALSE), ]
head(sig2)

pdf("~/repo/eQTL-2D/docs/manuscript/images/3dplots.pdf", width=15, height=15)
grid.arrange(
	plot3D(sig2[1,], xmat, resphen)[[3]],
	plot3D(sig2[2,], xmat, resphen)[[3]],
	plot3D(sig2[3,], xmat, resphen)[[3]],
	plot3D(sig2[4,], xmat, resphen)[[3]],
	plot3D(sig2[5,], xmat, resphen)[[3]],
	plot3D(sig2[6,], xmat, resphen)[[3]],
	plot3D(sig2[7,], xmat, resphen)[[3]],
	plot3D(sig2[8,], xmat, resphen)[[3]],
	plot3D(sig2[9,], xmat, resphen)[[3]],
	ncol=3
)
dev.off()

####################
# Interaction type #
####################

xtable(as.data.frame(a <- table(sig$vc)))
mat <- rbind(a, rep(mean(a),4))
chisq.test(mat, simulate.p.value = TRUE, B = 10000)
with(sig, tapply(pint, list(vc, type), mean))

xtable(as.data.frame(table(sig$type)))


##########################################
# Proportion of additive vs non-additive #
##########################################


prop <- subset(sig, select=c(varA, varI))
prop <- prop[order(prop$varA+prop$varI, decreasing=T),]
prop$index <- 1:nrow(prop)
prop <- melt(prop, id=c("index"))
prop <- prop[nrow(prop):1, ]
prop$variable <- factor(prop$variable, levels=c("varI", "varA"))
levels(prop$variable) <- c("Non-additive", "Additive")
ggplot(prop, aes(y=value, x=index)) + 
	geom_bar(stat="identity", width=1, size=0.1, colour="black", aes(fill=variable)) +
	scale_fill_brewer("Proportion of phenotypic variance", type="qual") +
	ylab("Phenotypic variance") +
	coord_flip()
ggsave("~/repo/eQTL-2D/docs/manuscript/images/proportion_additive.pdf", width=10, height=10)


##########################################
# Interaction effect vs allele frequency #
##########################################

maf <- apply(xmat, 2, function(x) {
	sum(x, na.rm=T) / (2*length(x))
})
sig$maf1 <- maf[sig$pos1]
sig$maf2 <- maf[sig$pos2]
sig$mmaf <- (sig$maf1 + sig$maf2) / 2
head(sig)
with(sig, tapply(maf2, vc, mean))
# There is no effect of MAF on which VC is most significant



###############
# Replication #
###############

load("~/repo/eQTL-2D/replication/run/gg_replication.RData")
a <- qqPlot(10^-newsig$replication_pnest)

pdf(file="~/repo/eQTL-2D/docs/manuscript/images/replication_qqplot.pdf")
plot(a, xlim=c(0, max(a)), ylim=c(0, max(a)), main="Replication p-values")
abline(a=0, b=1)
dev.off()




# Plot the Bonferroni significant probes

sig$probeid <- match(sig$probename, colnames(resphen))

a <- subset(bsgs, probegene == "TMEM149")

l <- list()
for(i in 1:nrow(a))
{
	l[[i]] <- plot3D(a[i,], xmat, resphen)[[3]]
}


b <- subset(bsgs, probegene == "MBNL1")

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



sig2 <- subset(sig, rep == "Bonf")
sig2$probeid <- match(sig2$probename, colnames(resphen))
grid.arrange(
	plot3D(sig2[4,], xmat, resphen)[[3]],
	plot3D(sig2[5,], xmat, resphen)[[3]],
	plot3D(sig2[6,], xmat, resphen)[[3]],
	plot3D(sig2[7,], xmat, resphen)[[3]],
	plot3D(sig2[8,], xmat, resphen)[[3]],
	plot3D(sig2[9,], xmat, resphen)[[3]],
	ncol=3
)



grid.arrange(
	plot3D(sig2[10,], xmat, resphen)[[3]],
	plot3D(sig2[11,], xmat, resphen)[[3]],
	plot3D(sig2[12,], xmat, resphen)[[3]],
	plot3D(sig2[13,], xmat, resphen)[[3]],
	plot3D(sig2[14,], xmat, resphen)[[3]],
	plot3D(sig2[15,], xmat, resphen)[[3]],
	plot3D(sig2[16,], xmat, resphen)[[3]],
	plot3D(sig2[17,], xmat, resphen)[[3]],
	ncol=3
)
