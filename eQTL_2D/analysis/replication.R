library(lattice)
library(latticeExtra)
library(ggplot2)

load("/Users/explodecomputer/git/wrayvisscher/eQTL_2D/analysis/residuals.RData")
load("/Users/explodecomputer/git/wrayvisscher/eQTL_2D/analysis/allsub.RData")
load("/Users/explodecomputer/git/wrayvisscher/eQTL_2D/analysis/clean_geno_final.RData")
load("/Users/explodecomputer/git/wrayvisscher/eQTL_2D/analysis/ggdata.RData")

gen[gen == "NC"] <- NA
gen[gen == "AA"] <- "0"
gen[gen == "AB"] <- "1"
gen[gen == "BB"] <- "2"
gen <- matrix(as.numeric(gen), nrow(gen), ncol(gen))

ls()


sig <- subset(allsub, pfull > 16.5 & propG > 0.05)
dups <- with(sig, paste(probe, chr1, chr2))
sig <- sig[!duplicated(dups), ]
sig <- subset(sig, snp1 != "rs11036212" & snp2 != "rs11036212" & probe != "ILMN_1688753")
dim(sig)

temp <- with(probeinfo, data.frame(probe=PROBE_ID, probechr=CHR))
sig <- merge(sig, temp, by="probe")
dim(sig)
head(sig)

dim(ciscis <- subset(sig, chr1 == probechr & chr2 == probechr))
dim(cistrans <- subset(sig, (chr1 == probechr & chr2 != probechr) | (chr1 != probechr & chr2 == probechr)))
dim(transtrans <- subset(sig, chr1 != probechr & chr2 != probechr))

sig$type <- NA
sig$type[with(sig, chr1 == probechr & chr2 == probechr)] <- "cis-cis"
sig$type[with(sig, (chr1 == probechr & chr2 != probechr) | chr1 != probechr & chr2 == probechr)] <- "cis-trans"
sig$type[with(sig, chr1 != probechr & chr2 != probechr)] <- "trans-trans"

head(sig)

lapply(sig, class)
sig$snp1 <- as.character(sig$snp1)
sig$snp2 <- as.character(sig$snp2)

allsnps <- with(sig, unique(c(snp1, snp2)))

length(allsnps)

table(allsnps %in% snp$Name)


table(sig$snp1 %in% snp$Name, sig$snp2 %in% snp$Name)


sig2 <- subset(sig, sig$snp1 %in% snp$Name & sig$snp2 %in% snp$Name)
dim(sig2)


head(sig2)




replicate <- function(
	prinfo,
	probe,
	gen,
	snp,
	row
) {
	# Get the row name in probe
	# Get the two SNP rows in gen

	prow <- which(prinfo$PROBE_ID == row$probe[1])
	stopifnot(length(prow) == 1)

	srow1 <- which(snp$Name == row$snp1)
	srow2 <- which(snp$Name == row$snp2)
	stopifnot(length(srow1) == 1)
	stopifnot(length(srow2) == 1)

	a <- anova(lm(probe[prow, ] ~ as.factor(gen[srow1, ])*as.factor(gen[srow2, ])))$P[1]
	return(a)
}


replicate(prinfo, probe, gen, snp, sig2[4,])

rep <- array(0, nrow(sig2))
for(i in 1:nrow(sig2))
{
	rep[i] <- replicate(prinfo, probe, gen, snp, sig2[i,])

}


# Of the 52 remaining, 3 interactions replicate
# 2 are cis-trans, one is cis-cis
sig2[which(rep * length(rep) < 0.05), ]
rep[which(rep * length(rep) < 0.05)]

sigrep <- sig2[which(rep * length(rep) < 0.05), ]


scrutinise <- function(res, geno, phen, z=45)
{
	a <- tapply(phen[,res$probeid], list(geno[,res$pos1], geno[, res$pos2]), mean)
	b <- table(geno[,res$pos1], geno[, res$pos2])
	p <- cloud(a, panel.3d.cloud=panel.3dbars, col="black", col.facet=c("#e5f5e0", "#A1D99B", "#31A354"), xbase=0.9, ybase=0.9,
	xlab="SNP1", ylab="SNP2", zlab="Phenotype",
	default.scales=list(arrows=F),
	screen = list(z = z, x = -60, y = 3))

	return(list(a,b,p))
}


scrutinise_rep <- function(
	prinfo,
	probe,
	gen,
	snp,
	row,
	z=45
) {
	# Get the row name in probe
	# Get the two SNP rows in gen

	prow <- which(prinfo$PROBE_ID == row$probe[1])
	stopifnot(length(prow) == 1)

	srow1 <- which(snp$Name == row$snp1)
	srow2 <- which(snp$Name == row$snp2)
	stopifnot(length(srow1) == 1)
	stopifnot(length(srow2) == 1)


	a <- tapply(probe[prow,], list(gen[srow1, ], gen[srow2, ]), mean)
	b <- table(gen[srow1, ], gen[srow2, ])
	p <- cloud(a, panel.3d.cloud=panel.3dbars, col="black", col.facet=c("#e5f5e0", "#A1D99B", "#31A354"), xbase=0.9, ybase=0.9,
	xlab="SNP1", ylab="SNP2", zlab="Phenotype",
	default.scales=list(arrows=F),
	screen = list(z = z, x = -60, y = 3))

	return(list(a,b,p))
}


scrutinise(sigrep[1,], xmat, resphen)
scrutinise(sigrep[3,], xmat, resphen, z=-30)

dev.new()
scrutinise_rep(prinfo, probe, gen, snp, sigrep[1,])
dev.new()
scrutinise_rep(prinfo, probe, gen, snp, sigrep[3,])


sigrep

original pval
rep pval
gene
chr
snp1
maf1
snp2
maf2

snp$maf <- apply(gen, 1, function(x) {
	sum(x, na.rm=T) / (2*sum(!is.na(x)))
	})
index <- snp$maf > 0.5
snp$maf[index] <- 1 - snp$maf[index]
hist(snp$maf)
table(snp$maf > 0.05)

bim$maf <- apply(xmat, 2, function(x) {
	sum(x, na.rm=T) / (2*sum(!is.na(x)))
	})
index <- bim$maf > 0.5
bim$maf[index] <- 1 - bim$maf[index]
hist(bim$maf)
table(bim$maf > 0.05)

probeinfo[sigrep$probeid,]


subset(snp, Name %in% c("rs7985085", "rs2241623"))


cor(xmat[, 488075], xmat[, 488068])


subset(snp, Name %in% c("rs10847601", "rs229670"))
subset(bim, V2 %in% c("rs10847601", "rs229670"))







