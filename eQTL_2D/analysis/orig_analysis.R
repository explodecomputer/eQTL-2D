library(lattice)
library(latticeExtra)

cloud(tab, panel.3d.cloud=panel.3dbars, col="black", col.facet=2:4, xbase=0.9, ybase=0.9)

sig <- subset(allsub, pfull > 16.5 & propG > 0.05)
dups <- with(sig, paste(probe, chr1, chr2))
sig <- sig[!duplicated(dups), ]

load("/Users/ghemani/wrayvisscher/eQTL_2D/data/residuals.RData")
load("/Users/ghemani/wrayvisscher/eQTL_2D/v4/allsub.RData")
load("/Users/ghemani/wrayvisscher/eQTL_2D/data/clean_geno_final.RData")

temp <- with(probeinfo, data.frame(probe=PROBE_ID, probechr=CHR))
sig <- merge(sig, temp, by="probe")

dim(sig)
dim(ciscis <- subset(sig, chr1 == probechr & chr2 == probechr))
dim(cistrans <- subset(sig, (chr1 == probechr & chr2 != probechr) | (chr1 != probechr & chr2 == probechr)))
dim(transtrans <- subset(sig, chr1 != probechr & chr2 != probechr))

sig$type <- NA
sig$type[with(sig, chr1 == probechr & chr2 == probechr)] <- "cis-cis"
sig$type[with(sig, (chr1 == probechr & chr2 != probechr) | chr1 != probechr & chr2 == probechr)] <- "cis-trans"
sig$type[with(sig, chr1 != probechr & chr2 != probechr)] <- "trans-trans"


# 140 significant
# 71 probes with significant hits

# cis-cis: 19
# cis-trans: 95
# trans-trans: 24


with(sig, tapply(pfull, type, mean))
library(ggplot2)

ggplot(sig, aes(x=propG)) + geom_histogram() + facet_grid(.~type)

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

scrutinise(sig[4,], xmat, resphen, z=130)
scrutinise(sig[5,], xmat, resphen)
scrutinise(sig[6,], xmat, resphen, z=20)
scrutinise(sig[7,], xmat, resphen)
scrutinise(sig[8,], xmat, resphen)
scrutinise(sig[9,], xmat, resphen)


p1 <- scrutinise(sig[1,], xmat, resphen) #
p2 <- scrutinise(sig[2,], xmat, resphen) #
p3 <- scrutinise(sig[3,], xmat, resphen, z=30) #
p4 <- scrutinise(sig[10,], xmat, resphen) #

pdf("epistasis_examples.pdf", width=20, height=20)
print(p1[[3]], position=c(0, 0.5, 0.5, 1), more=T)
print(p2[[3]], position=c(0.5, 0.5, 1, 1), more=T)
print(p3[[3]], position=c(0, 0, 0.5, 0.5), more=T)
print(p4[[3]], position=c(0.5, 0, 1, 0.5), more=T)
dev.off()

p1[[2]]
p2[[2]]
p3[[2]]
p4[[2]]

