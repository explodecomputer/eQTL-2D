library(reshape2)
library(ggplot2)


#================================================================================#
#================================================================================#


summariseVg <- function(dat)
{
	dat <- as.data.frame(dat)
	dat$varA <- dat$a. + dat$.a
	dat$varD <- dat$d. + dat$.d
	dat$varI <- dat$aa + dat$ad + dat$da + dat$dd
	return(dat)	
}

standardiseVg <- function(dat)
{
	vG <- apply(dat, 1, sum)
	dat <- dat / vG
	dat <- summariseVg(dat)
	return(dat)
}


sortVar <- function(bsgs, fehr, egcut, cohort)
{
	varA <- data.frame(bsgs = bsgs$varA, fehr = fehr$varA, egcut = egcut$varA)
	varD <- data.frame(bsgs = bsgs$varD, fehr = fehr$varD, egcut = egcut$varD)
	varI <- data.frame(bsgs = bsgs$varI, fehr = fehr$varI, egcut = egcut$varI)

	varA$index <- 1:nrow(varA)
	varA <- varA[order(varA[[cohort]], decreasing=T), ]
	varD <- varD[varA$index, ]
	varI <- varI[varA$index, ]
	varA <- subset(varA, select=-c(index))
	vA <- melt(varA)
	vA$vc <- "A"
	vD <- melt(varD)
	vD$vc <- "D"
	vI <- melt(varI)
	vI$vc <- "I"

	vA$index <- rep(1:nrow(bsgs), 3)
	vD$index <- rep(1:nrow(bsgs), 3)
	vI$index <- rep(1:nrow(bsgs), 3)

	dat <- rbind(vA, vD, vI)
	names(dat) <- c("dataset", "Variance", "Component", "index")
	levels(dat$dataset) <- c("BSGS", "Fehrmann", "EGCUT")

	return(dat)
}


sortVar2 <- function(bsgs, fehr, egcut)
{
	bsgs <- bsgs[order(bsgs$varA, decreasing=FALSE), ]
	fehr <- fehr[order(fehr$varA, decreasing=FALSE), ]
	egcut <- egcut[order(egcut$varA, decreasing=FALSE), ]	
	
	varA <- data.frame(bsgs = bsgs$varA, fehr = fehr$varA, egcut = egcut$varA)
	varD <- data.frame(bsgs = bsgs$varD, fehr = fehr$varD, egcut = egcut$varD)
	varI <- data.frame(bsgs = bsgs$varI, fehr = fehr$varI, egcut = egcut$varI)

	vA <- melt(varA)
	vA$vc <- "A"
	vD <- melt(varD)
	vD$vc <- "D"
	vI <- melt(varI)
	vI$vc <- "I"

	vA$index <- rep(1:nrow(bsgs), 3)
	vD$index <- rep(1:nrow(bsgs), 3)
	vI$index <- rep(1:nrow(bsgs), 3)

	dat <- rbind(vA, vD, vI)
	names(dat) <- c("dataset", "Variance", "Component", "index")
	levels(dat$dataset) <- c("BSGS", "Fehrmann", "EGCUT")

	return(dat)
}


plotvarG <- function(dat)
{
	dat$Component <- factor(dat$Component, levels=c("D", "I", "A"))
	levels(dat$Component) <- c("Dominance", "Interaction", "Additive")
	p1 <- ggplot(dat, aes(y=Variance, x=index)) +
		geom_bar(stat="identity", aes(fill=Component, colour=Component), position=position_stack(width=0)) +
		scale_fill_brewer("Variance component") +
		scale_colour_brewer("Variance component") +
		ylab("Genetic variance") + 
		xlab("") +
		facet_grid(dataset ~ .)
	return(p1)
}


#================================================================================#
#================================================================================#


load("~/repo/eQTL-2D/analysis/interaction_list_replication_summary.RData")


#================================================================================#
#================================================================================#


sig$type <- "cis-cis"
sig$type[with(sig, chr1 == probechr & chr2 != probechr)] <- "cis-trans"
sig$type[with(sig, chr1 != probechr & chr2 == probechr)] <- "cis-trans"
sig$type[with(sig, chr1 != probechr & chr2 != probechr)] <- "trans-trans"

bsgs <- do.call(rbind, sig$vc)
fehr <- do.call(rbind, sig$vc_fehr)
egcut <- do.call(rbind, sig$vc_egcut)

bsgs <- summariseVg(bsgs)
fehr <- summariseVg(fehr)
egcut <- summariseVg(egcut)

dat <- sortVar(bsgs, fehr, egcut, "fehr")
dat2 <- sortVar2(bsgs, fehr, egcut)

dat3 <- dat2[with(dat2, order(dataset, Component, Variance)), ]
dat3$index <- rep(1:501, 9)


#================================================================================#
#================================================================================#


plotvarG(dat2)
ggsave(file="~/repo/eQTL-2D/analysis/images/compare_vc.pdf", width=7, height=7)

ggplot(dat3, aes(y=Variance, x=index)) +
	geom_bar(stat="identity", position=position_stack(width=0)) +
	ylab("Proportion of phenotypic variance explained") + 
	xlab("") +
	facet_grid(Component ~ dataset, scales="free_y")
ggsave(file="~/repo/eQTL-2D/analysis/images/compare_vc2.pdf", width=7, height=7)


#================================================================================#
#================================================================================#


ggplot(dat2, aes(x=Variance)) + geom_histogram(binwidth=0.01) + facet_grid(dataset ~ Component)



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



# Proportion of SNPs with known marginal effects
with(sig, table(is.na(marginal_gene1), is.na(marginal_gene2)))




biggestVC <- function(vc)
{
	a <- apply(subset(vc, select=c(aa, ad, da, dd)), 1, which.max)
	a[a == 1] <- "aa"
	a[a == 2] <- "ad"
	a[a == 3] <- "ad"
	a[a == 4] <- "dd"
	vc$mainI <- a
	return(vc)
}

bsgs <- biggestVC(bsgs)

tab <- table(bsgs$mainI)
chisq.test(x=tab / sum(tab), y=c(0.25, 0.5, 0.25))



bsgs$varG <- with(bsgs, sum(varA, varD, varI))

s <- apply(subset(bsgs, select=c(varA, varD, varI)), 2, sum)
s / sum(s)


