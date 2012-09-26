library(lattice)
library(latticeExtra)
library(ggplot2)

load("/Users/ghemani/wrayvisscher/eQTL_2D/data/residuals.RData")
load("/Users/ghemani/wrayvisscher/eQTL_2D/v4/allsub.RData")
load("/Users/ghemani/wrayvisscher/eQTL_2D/data/clean_geno_final.RData")
load("/Users/ghemani/wrayvisscher/eQTL_2D/data/ggdata.RData")

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




}










