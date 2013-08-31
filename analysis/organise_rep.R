library(noia)
library(plyr)

##


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

calcGenoFreq <- function(geno)
{
	dat <- data.frame(snp = colnames(geno), freq = apply(geno, 2, function(x)
	{
		sum(x, na.rm=T) / (2 * length(na.omit(x)))
	}))
	return(dat)
}

calcGenoFreqAll <- function(genolist)
{
	l <- genolist
	n <- length(l)

	dat <- calcGenoFreq(l[[1]])

	for(i in 2:n)
	{
		dat <- merge(dat, calcGenoFreq(l[[i]]), by="snp")
	}
	names(dat) <- c("SNP", paste("f", 1:n, sep=""))
	return(dat)
}

doTest <- function(geno, phen, bsig, row, thresh = 0.01)
{
	snp1 <- geno[, which(colnames(geno) == bsig$snp1[row])]
	snp2 <- geno[, which(colnames(geno) == bsig$snp2[row])]
	y <- phen[, which(colnames(phen) == bsig$probename[row])]

	mod <- linearRegression(y, cbind(snp1, snp2) + 1)


	index <- c(2,4,3,7,5,8,6,9)

	# get effect and variance and pval

	dat <- data.frame(vc = names(mod$E)[index], E = mod$E[index], V = mod$variances[index], p = mod$pvalues[index], sig = mod$pvalues[index] < thresh)
	return(dat)
}

compareTests <- function(genolist, phenlist, bsig, row)
{
	l <- list()
	l[[1]] <- doTest(genolist[[1]], phenlist[[1]], bsig, row)
	l[[2]] <- doTest(genolist[[2]], phenlist[[2]], bsig, row)
	l[[3]] <- doTest(genolist[[3]], phenlist[[3]], bsig, row)

	a <- matrix(NA, 2, 3)
	for(i in 1:3)
	{
		m <- l[[i]]

		a[, i] <- sign(m$E[1:2])
		a[, i] <- a[, i] * m$sig[1:2]
		
	}
	return(a)
}

Mode <- function(x) {
	ux <- unique(x)
	ux[which.max(tabulate(match(x, ux)))]
}

findSnpToFlip <- function(l, s)
{
	# SNP 1
	i <- l[s, ] == 0
	f <- 0

	if(sum(i) == 0)
	{
		m <- Mode(l[s, ])
		if(! all(l[s, ] == m) )
		{

			f <- which(l[s, ] != m)
		}
	}

	if(sum(i) == 1 & 1 %in% l[s, ] & -1 %in% l[s, ])
	{
		f <- which(l[s, ] == -1)
	}

	return(f)

}

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

makeFlipList <- function(genolist, phenlist, bsig)
{
	l <- list()
	for(i in 1:30)
	{
		l[[i]] <- compareTests(genolist, phenlist, bsig, i)
	}
	return(l)
}

doTheFlipping <- function(genolist, phenlist, bsig)
{
	l <- list()
	for(i in 1:nrow(bsig))
	{
		l[[i]] <- compareTests(genolist, phenlist, bsig, i)
		for(j in 1:2)
		{
			f <- findSnpToFlip(l[[i]], j)
			if(f != 0)
			{
				print(c(i,j))
				print(l[[i]])
				genolist[[f]] <- flip(bsig, i, paste("snp", j, sep=""), genolist[[f]])
			}
		}
	}
	return(genolist)
}

compareTestsEpi <- function(genolist, phenlist, bsig, row, thresh)
{
	l <- list()
	l[[1]] <- doTest(genolist[[1]], phenlist[[1]], bsig, row, thresh)
	l[[2]] <- doTest(genolist[[2]], phenlist[[2]], bsig, row, thresh)
	l[[3]] <- doTest(genolist[[3]], phenlist[[3]], bsig, row, thresh)

	a <- matrix(NA, 4, 3)

	for(i in 1:3)
	{
		m <- l[[i]]
		a[, i] <- sign(m$E[5:8]) * m$sig[5:8]
	}
	return(a)
}

compareTestsEpi4 <- function(genolist, phenlist, bsig, row, thresh)
{
	l <- list()
	l[[1]] <- doTest(genolist[[1]], phenlist[[1]], bsig, row, thresh)
	l[[2]] <- doTest(genolist[[2]], phenlist[[2]], bsig, row, thresh)
	l[[3]] <- doTest(genolist[[3]], phenlist[[3]], bsig, row, thresh)

	a <- matrix(NA, 4, 3)

	for(i in 1:3)
	{
		m <- l[[i]]
		a[, i] <- sign(m$E[5:8]) * m$sig[5:8]
	}

	b <- c(
		sum(a[,1] == a[,2]),
		sum(a[,1] == a[,3])
	)

	return(b)
}

compareTestsEpiLargest <- function(genolist, phenlist, bsig, row, thresh)
{
	l <- list()
	l[[1]] <- doTest(genolist[[1]], phenlist[[1]], bsig, row, thresh)
	l[[2]] <- doTest(genolist[[2]], phenlist[[2]], bsig, row, thresh)
	l[[3]] <- doTest(genolist[[3]], phenlist[[3]], bsig, row, thresh)

	a <- array(NA, 3)

	j <- which.min(l[[1]]$p[5:8])
	a[1] <- sign(l[[1]]$E[j]) * l[[1]]$sig[j]

	for(i in 2:3)
	{
		m <- l[[i]]
		a[i] <- sign(m$E[j]) * m$sig[j]
	}
	return(a)
}

doAllTests <- function(genolist, phenlist, bsig, thresh)
{
	l <- list()
	for(i in 1:nrow(bsig))
	{
		l[[i]] <- list()
		l[[i]][[1]] <- doTest(genolist[[1]], phenlist[[1]], bsig, i, thresh)
		l[[i]][[2]] <- doTest(genolist[[2]], phenlist[[2]], bsig, i, thresh)
		l[[i]][[3]] <- doTest(genolist[[3]], phenlist[[3]], bsig, i, thresh)
	}
	return(l)
}

makeNoiaD <- function(noia)
{
	noiad <- list()
	k <- 1
	for(i in 1:length(noia))
	{
		for(j in 1:3)
		{
			noiad[[k]] <- data.frame(noia[[i]][[j]], interaction = i, dataset = j)
			k <- k + 1
		}
	}
	noiad <- rbind.fill(noiad)
	return(noiad)
}

compareTestsEpiLargest <- function(genolist, phenlist, bsig, row, thresh)
{
	l <- list()
	l[[1]] <- doTest(genolist[[1]], phenlist[[1]], bsig, row, thresh)
	l[[2]] <- doTest(genolist[[2]], phenlist[[2]], bsig, row, thresh)
	l[[3]] <- doTest(genolist[[3]], phenlist[[3]], bsig, row, thresh)

	a <- array(NA, 3)

	j <- which.min(l[[1]]$p[5:8])
	a[1] <- sign(l[[1]]$E[j]) * l[[1]]$sig[j]

	for(i in 2:3)
	{
		m <- l[[i]]
		a[i] <- sign(m$E[j]) * m$sig[j]
	}
	return(a)
}

compareTestsEpiLargestAll <- function(genolist, phenlist, bsig, thresh)
{
	a <- matrix(NA, nrow(bsig), 3)
	for(i in 1:nrow(bsig))
	{
		a[i,] <- compareTestsEpiLargest(genolist, phenlist, bsig, i, thresh)
	}
	return(a)
}

compareTestsEpi4All <- function(genolist, phenlist, bsig, thresh)
{
	a <- matrix(NA, nrow(bsig), 2)
	for(i in 1:nrow(bsig))
	{
		a[i,] <- compareTestsEpi4(genolist, phenlist, bsig, i, thresh)
	}
	return(a)
}

compareTestsEpiLargest2 <- function(genolist, phenlist, bsig, row, thresh)
{
	l <- list()
	l[[1]] <- doTest(genolist[[1]], phenlist[[1]], bsig, row, thresh)
	l[[2]] <- doTest(genolist[[2]], phenlist[[2]], bsig, row, thresh)
	l[[3]] <- doTest(genolist[[3]], phenlist[[3]], bsig, row, thresh)

	a <- array(NA, 3)

	j <- which.min(l[[1]]$p[5:8])
	a <- sign(l[[1]]$E[j])

	b <- array(NA, 2)

	for(i in 2:3)
	{
		m <- l[[i]]
		b[i-1] <- sign(m$E[j]) == a & which.min(m$p[5:8])[1] == j
	}
	print(b)
	return(b)
}

compareTestsEpiLargest2All <- function(genolist, phenlist, bsig, thresh)
{
	a <- matrix(NA, nrow(bsig), 2)
	for(i in 1:nrow(bsig))
	{
		a[i,] <- compareTestsEpiLargest2(genolist, phenlist, bsig, i, thresh)
	}
	return(a)
}

compareTestsEpiLargest1 <- function(genolist, phenlist, bsig, row, thresh)
{
	l <- list()
	l[[1]] <- doTest(genolist[[1]], phenlist[[1]], bsig, row, thresh)
	l[[2]] <- doTest(genolist[[2]], phenlist[[2]], bsig, row, thresh)
	l[[3]] <- doTest(genolist[[3]], phenlist[[3]], bsig, row, thresh)

	a <- array(NA, 3)

	j <- which.min(l[[1]]$p[5:8])
	a <- sign(l[[1]]$E[j])

	b <- array(NA, 2)

	for(i in 2:3)
	{
		m <- l[[i]]
		b[i-1] <- sign(m$E[j]) == a
	}
	print(b)
	return(b)
}

compareTestsEpiLargest1All <- function(genolist, phenlist, bsig, thresh)
{
	a <- matrix(NA, nrow(bsig), 2)
	for(i in 1:nrow(bsig))
	{
		a[i,] <- compareTestsEpiLargest1(genolist, phenlist, bsig, i, thresh)
	}
	return(a)
}


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

b_phen <- getPhen(b_mod, b_newsig)
e_phen <- getPhen(e_mod, e_newsig)
f_phen <- getPhen(f_mod, f_newsig)
phenlist <- list(b_phen, e_phen, f_phen)

b_geno <- getGeno(b_mod, b_newsig)
e_geno <- getGeno(e_mod, e_newsig)
f_geno <- getGeno(f_mod, f_newsig)
genolist <- list(b_geno, e_geno, f_geno)

f <- calcGenoFreqAll(genolist)
cor(f[,2:4])


##

load("~/repo/eQTL-2D/analysis/interaction_list_meta_analysis.RData")
sig <- subset(meta, filter !=3)
bsig <- subset(sig, pnest_meta > -log10(0.05/434))
bsig <- bsig[order(bsig$probegene), ]


##


a <- makeFlipList(genolist, phenlist, bsig)
a

gl2 <- doTheFlipping(genolist, phenlist, bsig)

gl3 <- doTheFlipping(genolist, phenlist, sig)


f1 <- calcGenoFreqAll(genolist)
pairs(f[,2:4])

f2 <- calcGenoFreqAll(gl2)
dev.new()
pairs(f2[,2:4])

f3 <- calcGenoFreqAll(gl3)
dev.new()
pairs(f3[,2:4])


##


# e_geno2 <- e_geno
# e_geno <- flip(bsig, 2, "snp2", e_geno)
# e_geno <- flip(bsig, 3, "snp2", e_geno)
# e_geno <- flip(bsig, 6, "snp1", e_geno)
# e_geno <- flip(bsig, 11, "snp2", e_geno)
# e_geno <- flip(bsig, 12, "snp2", e_geno)
# e_geno <- flip(bsig, 13, "snp2", e_geno)
# e_geno <- flip(bsig, 15, "snp2", e_geno)
# e_geno <- flip(bsig, 22, "snp2", e_geno)

# b_geno2 <- b_geno
# b_geno <- flip(bsig, 5, "snp2", b_geno)
# b_geno <- flip(bsig, 22, "snp2", b_geno)

# f_geno2 <- f_geno
# f_geno <- flip(bsig, 10, "snp2", f_geno)
# f_geno <- flip(bsig, 22, "snp2", f_geno)



##

genolist <- gl3
save(genolist, phenlist, file="bsgs_egcut_fehr_data.RData")

noia <- doAllTests(genolist, phenlist, sig, 1)
bnoia <- doAllTests(genolist, phenlist, bsig, 1)
noiad <- makeNoiaD(noia)
bnoiad <- makeNoiaD(bnoia)
save(noiad, bnoiad, file = "~/repo/eQTL-2D/analysis/noiad.RData")

