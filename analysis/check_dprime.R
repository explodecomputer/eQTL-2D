load("~/repo/eQTL-2D/analysis/interaction_list_meta_analysis.RData")
load("~/repo/eQTL-2D/data/clean_geno_final.RData")


library(snpStats)
genoToSnpMatrix <- function(geno, fam, bim)
{
	geno <- geno + 1
	a <- matrix(as.raw(geno), nrow(geno), ncol(geno))
	rownames(a) <- fam$iid
	colnames(a) <- bim$snpname
	snpmatrix <- new("SnpMatrix", .Data=a)
	return(snpmatrix)
}

#============================================================#
#============================================================#


dim(meta)
sig <- subset(meta, filter != 3)
dim(sig)

names(sig)

sig$type <- "cis-cis"
sig$type[sig$chr1 != sig$probechr & sig$chr2 == sig$probechr] <- "cis-trans"
sig$type[sig$chr1 == sig$probechr & sig$chr2 != sig$probechr] <- "cis-trans"
sig$type[sig$chr1 != sig$probechr & sig$chr2 != sig$probechr] <- "trans-trans"

table(sig$type)

ciscis <- subset(sig, type=="cis-cis")
ciscis$position2[36] <- 23299135
ciscis$position2[37] <- 59874129
ciscis$gap <- abs(ciscis$position1 - ciscis$position2) / 1000000

ciscis


#============================================================#
#============================================================#


index <- bim$V2 %in% with(ciscis, c(snp1, snp2))
table(index)

sb <- bim[index, ]
sg <- xmat[, index]

dim(sb)
dim(sg)

gr <- genoToSnpMatrix(sg, fam, sb)


dprime <- ld(gr, depth=nrow(sb), stats="D.prime")

class(dprime)
dprime
max(dprime)
table(dprime > 0.3)

l <- array(0, nrow(ciscis))
for(i in 1:nrow(ciscis))
{
	snp <- with(ciscis, c(snp1[i], snp2[i]))
	a <- gr[, which(sb$V2 %in% snp)]
	stopifnot(ncol(a)==2)
	l[i] <- ld(a, depth=1, stats="D.prime")[3]
	cat(i, " - ", l[i], "\n")
}

