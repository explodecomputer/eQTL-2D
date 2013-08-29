library(noia)
library(plyr)

##


load("~/repo/eQTL-2D/replication/results/replication_EGCUT.RData")
e_gcm <- gcm
e_gcs <- gcs
e_mod <- mod
e_newsig <- newsig

load("~/repo/eQTL-2D/replication/results/replication_GrngHT12v3.RData")
f_gcm <- gcm
f_gcs <- gcs
f_mod <- mod
f_newsig <- newsig

load("~/repo/eQTL-2D/replication/results/original_bsgs.RData")
b_gcm <- gcm
b_gcs <- gcs
b_mod <- mod
b_newsig <- newsig

rm(mod, gcm, gcs, newsig)


##


## get phenotypes (standardised)
## get genotypes (SNP names, correct allele frequency)


normalise <- function(x)
{
	(x - mean(x, na.rm=T)) / sd(x, na.rm=T)
}

getPhen <- function(mod, newsig)
{
	out <- llply(mod, function(x)
	{
		if(!is.na(x))
		{
			return(normalise(x$phen))
		} else {
			return(NA)
		}
	})
	out <- do.call(cbind, out)
	colnames(out) <- as.character(newsig$probename)
	index <- duplicated(colnames(out))
	out <- out[, !index]
	return(out)
}


b_phen <- getPhen(b_mod, b_newsig)
dim(b_phen)

e_phen <- getPhen(e_mod, e_newsig)
dim(e_phen)

f_phen <- getPhen(f_mod, f_newsig)
dim(f_phen)



getGeno <- function(mod, newsig)
{
	out <- llply(mod, function(x)
	{
		return(x$gen - 1)
	})
	out <- do.call(cbind, out)
	snp1 <- as.character(newsig$snp1)
	snp2 <- as.character(newsig$snp2)
	allsnp <- array("0", 2*length(snp1))
	allsnp[seq(1, length(allsnp), 2)] <- snp1
	allsnp[seq(2, length(allsnp), 2)] <- snp2
	colnames(out) <- as.character(allsnp)
	index <- duplicated(colnames(out))
	out <- out[, !index]
	return(out)
}



test_getGeno <- function(mod, newsig, geno, i)
{
	j <- which(colnames(b_geno) == b_newsig$snp1[i])
	print(cor(b_mod[[i]]$gen[,1], b_geno[,j], use="pair"))

	j <- which(colnames(b_geno) == b_newsig$snp2[i])
	print(cor(b_mod[[i]]$gen[,2], b_geno[,j], use="pair"))
}

test_getGeno(b_mod, b_newsig, b_geno, 1)


b_geno <- getGeno(b_mod, b_newsig)
e_geno <- getGeno(e_mod, e_newsig)
f_geno <- getGeno(f_mod, f_newsig)
dim(b_geno)
dim(e_geno)
dim(f_geno)



calcGenoFreq <- function(geno)
{
	dat <- data.frame(snp = colnames(geno), freq = apply(geno, 2, function(x)
	{
		sum(x, na.rm=T) / (2 * length(na.omit(x)))
	}))
	return(dat)
}

calcGenoFreqAll <- function(...)
{
	l <- list(...)
	n <- length(l)

	dat <- calcGenoFreq(l[[1]])

	for(i in 2:n)
	{
		dat <- merge(dat, calcGenoFreq(l[[i]]), by="snp")
	}
	names(dat) <- c("SNP", paste("f", 1:n, sep=""))
	return(dat)
}

f <- calcGenoFreqAll(b_geno, f_geno, e_geno)
head(f)

cor(f[,2:4])
pairs(f[,2:4])



##

load("~/repo/eQTL-2D/analysis/interaction_list_meta_analysis.RData")
sig <- subset(meta, filter !=3)
bsig <- subset(sig, pnest_meta > -log10(0.05/434))
bsig <- bsig[order(bsig$probegene), ]
dim(bsig)


bsig

doTest <- function(geno, phen, bsig, row)
{
	snp1 <- geno[, which(colnames(geno) == bsig$snp1[row])]
	snp2 <- geno[, which(colnames(geno) == bsig$snp2[row])]
	y <- phen[, which(colnames(phen) == bsig$probename[row])]

	mod <- linearRegression(y, cbind(snp1, snp2) + 1)


	index <- c(2,4,3,7,5,8,6,9)

	# get effect and variance and pval

	dat <- data.frame(vc = names(mod$E)[index], E = mod$E[index], V = mod$variances[index], p = mod$pvalues[index], sig = mod$pvalues[index] < 0.01)
	return(dat)
}


doTest(b_geno, b_phen, bsig, 1)
doTest(f_geno, f_phen, bsig, 1)
doTest(e_geno, e_phen, bsig, 1)

doTest(b_geno, b_phen, bsig, 10)
doTest(f_geno, f_phen, bsig, 10)
doTest(e_geno, e_phen, bsig, 10)

doTest(b_geno, b_phen, bsig, 30)
doTest(f_geno, f_phen, bsig, 30)
doTest(e_geno, e_phen, bsig, 30)


compareTests <- function(row)
{
	l <- list()
	l[[1]] <- doTest(b_geno, b_phen, bsig, row)
	l[[2]] <- doTest(f_geno, f_phen, bsig, row)
	l[[3]] <- doTest(e_geno, e_phen, bsig, row)

	a <- matrix(NA, 2, 3)
	for(i in 1:3)
	{
		m <- l[[i]]

		a[, i] <- sign(m$E[1:2])
		a[, i] <- a[, i] * m$sig[1:2]
		
	}
	return(a)
}


makeFlipList <- function()
{
	l <- list()
	for(i in 1:30)
	{
		l[[i]] <- compareTests(i)
	}
	return(l)
}


a <- makeFlipList()



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



flipGeno <- function(x)
{
	x[x == 0] <- 3
	x[x == 2] <- 0
	x[x == 3] <- 2
	return(x)
}


flip <- function(bsig, row, snp, geno)
{
	snpid <- bsig[[snp]][row]
	i <- which(colnames(geno) == snpid)
	a <- flipGeno(geno[,i])
	print(cor(a, geno[,i], use="pair"))
	geno[,i] <- a
	return(geno)
}


e_geno2 <- e_geno

e_geno3 <- flip(bsig, 2, "snp2", e_geno)











