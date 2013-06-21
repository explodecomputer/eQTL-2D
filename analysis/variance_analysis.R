library(reshape2)
library(ggplot2)

load("~/repo/eQTL-2D/analysis/interaction_list_replication_summary.RData")


sig$type <- "cis-cis"
sig$type[with(sig, chr1 == probechr & chr2 != probechr)] <- "cis-trans"
sig$type[with(sig, chr1 != probechr & chr2 == probechr)] <- "cis-trans"
sig$type[with(sig, chr1 != probechr & chr2 != probechr)] <- "trans-trans"

for(i in 1:nrow(sig))
{
	sig$vc_fehr[[i]] <- sig$vc_fehr[[i]][[1]]
	sig$vc_egcut[[i]] <- sig$vc_egcut[[i]][[1]]
}

bsgs <- do.call(rbind, sig$vc)
fehr <- do.call(rbind, sig$vc_fehr)
egcut <- do.call(rbind, sig$vc_egcut)

standardiseVg <- function(dat)
{
	vG <- apply(dat, 1, sum)
	dat <- dat / vG
	dat <- as.data.frame(dat)
	dat$varA <- dat$a. + dat$.a
	dat$varD <- dat$d. + dat$.d
	dat$varI <- dat$aa + dat$ad + dat$da + dat$dd
	return(dat)
}

bsgs <- standardiseVg(bsgs)
fehr <- standardiseVg(fehr)
egcut <- standardiseVg(egcut)

varA <- data.frame(bsgs = bsgs$varA, fehr = fehr$varA, egcut = egcut$varA)
varD <- data.frame(bsgs = bsgs$varD, fehr = fehr$varD, egcut = egcut$varD)
varI <- data.frame(bsgs = bsgs$varI, fehr = fehr$varI, egcut = egcut$varI)

pairs(varA)
pairs(varD)
pairs(varI)

vA <- melt(varA)
vA$vc <- "A"
vD <- melt(varD)
vD$vc <- "D"
vI <- melt(varI)
vI$vc <- "I"

vA$index <- 1:nrow(vA)
vD$index <- 1:nrow(vD)
vI$index <- 1:nrow(vI)

dat <- rbind(vA, vD, vI)
dim(dat)
names(dat) <- c("dataset", "Variance", "Component", "index")
head(dat)

with(dat, tapply(Variance, list(dataset, Component), function(x) mean(x, na.rm=T)))



plotvarG <- function(dat)
{
	p1 <- ggplot(dat, aes(y=Variance, x=index)) +
		geom_bar(stat="identity", aes(fill=Component, colour=Component), position=position_stack(width=0)) +
		scale_fill_brewer("Variance component") +
		scale_colour_brewer("Variance component") +
		ylab("Phenotypic variance") + xlab("") +
		coord_flip() + theme(legend.position = "none") +
		# theme(axis.text.y=element_text(size=0), axis.ticks.y=element_line(size=0)) +
		facet_grid(. ~ dataset)
	return(p1)
}

plotvarG(dat)


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


