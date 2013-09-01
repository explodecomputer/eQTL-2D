library(plyr)

replicationTests <- function(geno, phen, ids)
{
	snp1 <- geno[, which(colnames(geno) == ids[[2]])[1]]
	snp2 <- geno[, which(colnames(geno) == ids[[3]])[1]]
	probe <- phen[, which(colnames(phen) == ids[[1]])[1]]

	# Statistical tests
	fullmod <- lm(phen ~ as.factor(snp1) * as.factor(snp2))
	margmod <- lm(phen ~ as.factor(snp1) + as.factor(snp2))
	names(inttest <- anova(margmod, fullmod))
	a <- inttest$Pr[2]
	p <- data.frame(c(ids, a))
	names(p) <- c("phen", "snp1", "snp2", "p")
	return(p)
}

doTest <- function(geno, phen, sig, row)
{

	ids <- list(sig$probename[row], sig$snp1[row], sig$snp2[row])
	p <- replicationTests(geno, phen, ids)
	return(p)
}

runTest <- function(geno, phen, sig)
{
	p <- list()
	for(i in 1:nrow(sig))
	{
		cat(i, "\n")
		p[[i]] <- doTest(geno, phen, sig, i)
	}
	p <- rbind.fill(p)
	return(p)
}




load("~/repo/eQTL-2D/data/bsgs_egcut_fehr_data.RData")
load("~/repo/eQTL-2D/data/bsgs_egcut_fehr_data.RData")
load("~/repo/eQTL-2D/analysis/interaction_list_meta_analysis.RData")
sig <- subset(meta, filter !=3)
bsig <- subset(sig, pnest_meta > -log10(0.05/434))
bsig <- bsig[order(bsig$probegene), ]

probes <- unique(as.character(sig$probename))
snps <- with(sig, unique(c(snp1, snp2)))

x1 <- match(snps, colnames(genolist[[1]]))
x2 <- match(snps, colnames(genolist[[2]]))
x3 <- match(snps, colnames(genolist[[3]]))

g1 <- genolist[[1]][,x1]
g2 <- genolist[[2]][,x2]
g3 <- genolist[[3]][,x3]
all(colnames(g1) == colnames(g2))
all(colnames(g1) == colnames(g3))

y1 <- match(probes, colnames(phenlist[[1]]))
y2 <- match(probes, colnames(phenlist[[2]]))
y3 <- match(probes, colnames(phenlist[[3]]))

p1 <- phenlist[[1]][,y1]
p2 <- phenlist[[2]][,y2]
p3 <- phenlist[[3]][,y3]
all(colnames(p1) == colnames(p2))
all(colnames(p1) == colnames(p3))

head(p1)


allg <- rbind(g1, g2, g3)
allp <- rbind(p1, p2, p3)
dim(allg)
dim(allp)

repg <- rbind(g2, g3)
repp <- rbind(p2, p3)
dim(repg)
dim(repp)


p_all <- runTest(allg, allp, sig)

p_bsgs <- runTest(genolist[[1]], phenlist[[1]], sig)

