library(plyr)

####

# Take a SNP / probe combination in replication dataset
# Choose another SNP randomly from remaining SNPs
# Check that it doesn't already exist as a real signal
# Perform test
# Do this 1000 times
# Keep the first 434


modifySig <- function(sig)
{
	sig$comb1 <- with(sig, paste(probename, snp1))
	sig$comb2 <- with(sig, paste(probename, snp2))
	sig$comb3 <- with(sig, paste(snp1, snp2))
	sig$comb4 <- with(sig, paste(snp2, snp1))
	return(sig)
}


chooseRandomSnpPair <- function(sigm)
{
	sig <- sigm
	flag <- 1
	while(flag == 1)
	{
		a <- sample(1:nrow(sig), 2, replace=FALSE)
		ran1 <- a[1]
		ran2 <- a[2]
		probename <- sig$probename[ran1]
		snp1 <- sig$snp1[ran1]
		snp2 <- sig$snp2[ran2]
		t1 <- paste(probename, snp1)
		t2 <- paste(probename, snp2)
		t3 <- paste(snp1, snp2)
		t4 <- paste(snp2, snp1)
		if(
			t2 %in% sig$comb2 |
			t3 %in% sig$comb3 |
			t4 %in% sig$comb4 |
			snp1 == snp2
		) {
			flag <- 1
		} else {
			flag <- 0
		}
	}
	return(list(probename, snp1, snp2))
}


replicationTests <- function(genolist, phenlist, ids)
{
	a <- list()

	for(i in 1:3)
	{
		geno <- genolist[[i]]
		probes <- phenlist[[i]]
		snp1 <- geno[, colnames(geno) == ids[[2]]]
		snp2 <- geno[, colnames(geno) == ids[[3]]]
		probe <- probes[, colnames(probes) == ids[[1]]]

		# Statistical tests
		fullmod <- lm(probe ~ as.factor(snp1) * as.factor(snp2))
		margmod <- lm(probe ~ as.factor(snp1) + as.factor(snp2))
		inttest <- anova(margmod, fullmod)
		a[[i]] <- inttest$P[2]
	}
	p <- data.frame(c(ids, a))
	names(p) <- c("probe", "snp1", "snp2", "bsgs", "egcut", "fehr")
	return(p)
}


doNullTest <- function(genolist, phenlist, sig)
{
	sigm <- modifySig(sig)
	ids <- chooseRandomSnpPair(sigm)
	p <- replicationTests(genolist, phenlist, ids)
	return(p)
}


runNull <- function(genolist, phenlist, sig, n)
{
	p <- list()
	for(i in 1:n)
	{
		p[[i]] <- doNullTest(genolist, phenlist, sig)
	}
	p <- rbind.fill(p)
	p$meta <- pchisq(-2*log(p$egcut) - 2*log(p$fehr), 4, lower.tail=FALSE)
	return(p)
}


##


load("~/repo/eQTL-2D/data/bsgs_egcut_fehr_data.RData")
load("~/repo/eQTL-2D/analysis/interaction_list_meta_analysis.RData")
sig <- subset(meta, filter !=3)
bsig <- subset(sig, pnest_meta > -log10(0.05/434))
bsig <- bsig[order(bsig$probegene), ]


p <- runNull(genolist, phenlist, sig, 434)
dev.new()
par(mfrow=c(2,2))
for(i in 1:4) hist(p[,i+3])
pnull <- p
save(pnull, file="~/repo/eQTL-2D/data/replication_null.RData")
