# Read in the epigpu output
# Calculate the variance components

datstat <- function(info, xmat, resphen)
{
	cat(" ... ",nrow(info), "rows to be analysed\n")
	a <- array(0, nrow(info))
	g <- array(0, nrow(info))
	for(i in 1:nrow(info))
	{
		if(i %% (nrow(info) / 100) == 0) cat(i/(nrow(info)/100)," ")
		l <- info[i,]
		if(is.na(l$chr1)) {
			g[i] <- NA
			a[i] <- NA
			next
		}
		mod <- anova(lm(resphen[,l$probeid] ~ xmat[,l$pos1] + xmat[,l$pos2] + as.factor(xmat[,l$pos1]) : as.factor(xmat[,l$pos2])))
		g[i] <- sum(mod$Sum[1:3]) / mod$Sum[4]
		a[i] <- sum(mod$Sum[1:2]) / sum(mod$Sum[1:3])

	}
	return(data.frame(g,a))
}

colCors = function(x, y) { 
	sqr = function(x) x*x
	if(!is.matrix(x)||!is.matrix(y)||any(dim(x)!=dim(y)))
	stop("Please supply two matrices of equal size.")
	x   = sweep(x, 2, colMeans(x, na.rm=T))
	y   = sweep(y, 2, colMeans(y, na.rm=T))
	cor = colSums(x*y, na.rm=T) /  sqrt(colSums(sqr(x), na.rm=T)*colSums(sqr(y), na.rm=T))
	return(cor)
}

calc.snpcor <- function(dat, geno) {
	snpnames <- colnames(geno)
	snp1pos <- match(dat$snp1, snpnames)
	snp2pos <- match(dat$snp2, snpnames)
	x <- geno[, snp1pos]
	y <- geno[, snp2pos]
	dat$snpcor <- colCors(x, y)
	return(dat)
}

full_vs_reduced <- function(dat, geno, phen) {
	dat$pnest <- NA
	for(i in 1:nrow(dat))
	{
		cat(i, "of", nrow(dat), "\n")
		x1 <- geno[, dat$pos1[i]]
		x2 <- geno[, dat$pos2[i]]
		y <- resphen[, dat$probeid[i]]
		full <- lm(y ~ as.factor(x1) + as.factor(x2) + as.factor(x1):as.factor(x2))
		reduced <- lm(y ~ as.factor(x1) + as.factor(x2))
		dat$pnest[i] <- -log10(anova(reduced, full)$P[2])
	}
	return(dat)
}

min.classsize <- function(dat, geno) {
	snpnames <- colnames(geno)
	snp1pos <- match(dat$snp1, snpnames)
	snp2pos <- match(dat$snp2, snpnames)
	x <- geno[, snp1pos]
	y <- geno[, snp2pos]
	xy <- x * 3 + y
	dat$minclasssize <- apply(xy, 2, function(x) min(table(x)))
	return(dat)
}

# /clusterdata/apps/R-2.14/bin/R --no-save --args results/result ../data/residuals.RData ../data/clean_geno_final.RData scratch/resphen ${start} ${end} 12 ../data/eqtl2d_objects.RData cols${i}.RData < collate_results.R


i <- as.numeric(commandArgs(T)[1])
rootres <- commandArgs(T)[2]
phenfile <- commandArgs(T)[3]
genofile <- commandArgs(T)[4]
threshold <- as.numeric(commandArgs(T)[5])
maxrsq <- as.numeric(commandArgs(T)[6])
minclass <- as.numeric(commandArgs(T)[7])
output <- commandArgs(T)[8]


# probe | probe number | hsq | complete | chr1 | chr2 | pos1 | pos2 | snp1 | snp2 | Pfull | Pint | df1 | df2 
# if probe has no values then single row with NA 


# Read in epiqpu output
load(paste(rootres, i, ".RData", sep=""))
dim(res)
if(nrow(res) == 1 & is.na(res$df1[1])) q()


# Remove non-autosomes
res <- subset(res, chr1 %in% 1:22 & chr2 %in% 1:22)
dim(res)
if(nrow(res) == 0) q()


# Filter based on 8df
res <- subset(res, df1 == 8)
dim(res)
if(nrow(res) == 0) q()


# Read in geno / pheno info
load(phenfile)
load(genofile)


# Filter based minimum number of individuals
res <- min.classsize(res, geno)
res <- subset(res, minclasssize >= minclass)
dim(res)
if(nrow(res) == 0) q()


# rsq < .1
res <- calc.snpcor(res, geno)
res <- subset(res, snpcor^2 <= maxrsq)
dim(res)
if(nrow(res) == 0) q()


# Choose the largest pfull for each chromosome pair
res <- res[order(res$chr1, res$chr2, res$pfull, decreasing=T), ]
res <- subset(res, !duplicated(paste(chr1, chr2)))


# Calculate additive vs non-additive
temp <- datstat(res, geno, resphen)
res$propG <- temp$g
res$propA <- temp$a


# Perform nested test
res <- full_vs_reduced(res, geno, phen)


save(res, file=paste(output, i, ".RData", sep=""))

