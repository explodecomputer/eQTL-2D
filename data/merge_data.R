mergeResiduals <- function()
{
	load("~/repo/eQTL-2D/data/residuals.RData")
	hsq1 <- hsq
	resphen1 <- resphen
	probeinfo1 <- probeinfo
	load("~/repo/eQTL-2D/data/residuals2.RData")
	hsq <- c(hsq1, hsq)
	resphen <- cbind(resphen1, resphen)
	probeinfo <- rbind(probeinfo1, probeinfo)

	save(hsq, probeinfo, resphen, file="~/repo/eQTL-2D/data/residuals_all.RData")
}

characteriseInteraction <- function(sig)
{
	sig$type <- NA
	index <- with(sig, chr1 == probechr & chr2 == probechr)
	sig$type[index] <- "cis-cis"
	index <- with(sig, (chr1 == probechr & chr2 != probechr) | (chr1 != probechr & chr2 == probechr))
	sig$type[index] <- "cis-trans"
	index <- with(sig, chr1 != probechr & chr2 != probechr)
	sig$type[index] <- "trans-trans"
	return(sig)
}

loadSigData <- function()
{
	load("~/repo/eQTL-2D/filtering/filtered_by_chr/sig_all_probes.RData")
	bim <- read.table("~/repo/eQTL-2D/data/clean_geno_final.bim", colClasses=c("character", "character", "numeric", "numeric", "character", "character"))
	sig$position1 <- bim$V4[sig$pos1]
	sig$position2 <- bim$V4[sig$pos2]
	sig$snp1 <- as.character(sig$snp1)
	sig$snp2 <- as.character(sig$snp2)
	sig$probename <- as.character(sig$probename)
	sig$probegene <- as.character(sig$probegene)

	# map probe chromosomes
	sig$probechr[sig$probegene == "MAP3K2"] <- 2
	sig$probechr[sig$probegene == "NAPSB"] <- 19
	sig$probechr[sig$probegene == "RPS26L"] <- 12
	sig$probechr[sig$probegene == "SAPS2"] <- 22
	sig$probechr[sig$probegene == "C15ORF28"] <- 22

	sig <- characteriseInteraction(sig)

	save(sig, file="~/repo/eQTL-2D/filtering/filtered_by_chr/sig_all_probes.RData")
}

temp <- mergeResiduals()
sig <- loadSigData()
