library(noia)
library(plyr)
library(biomaRt)

ReadOrig <- function(filename, objname, setname)
{
	load(filename)
	sig <- get(objname)

	sig <- subset(sig, select=c(chr1, chr2, snp1, snp2, pos1, pos2, probename, probegene, probechr, pfull, pnest, filter))
	sig$set <- setname
	sig$code <- with(sig, paste(probename, snp1, snp2))
	return(sig)
}

FormatRep <- function(sig)
{
	a <- sum(sig$replication_r^2 > 0.1)
	# cat(a, "pairs with high LD\n")

	b <- sum(sig$replication_nclass != 9)
	# cat(b, "pairs with < 9 classes\n")

	vc <- apply(subset(sig, replication_r^2 < 0.1 & replication_nclass == 9, select=c(.a, a., .d, d., aa, ad, da, dd)), 1, list)
	sig <- subset(sig, replication_r^2 < 0.1 & replication_nclass == 9, select=c(chr1, chr2, snp1, snp2, probename, probegene, probechr, replication_pfull, replication_pnest, filter, gcm, gcs))


	names(sig) <- c("chr1", "chr2", "snp1", "snp2", "probename", "probegene", "probechr", "pfull", "pnest", "filter", "gcm", "gcs")
	sig$vc <- vc
	sig$code <- with(sig, paste(probename, snp1, snp2))
	sig <- sig[order(sig$pnest, decreasing=T), ]
	return(sig)
}


posData <- function(sig)
{
	snplist <- unique(c(sig$snp1, sig$snp2))
	snpmart <- useMart("snp", dataset = "hsapiens_snp")
	snp_pos <- getBM(
		attributes = c("refsnp_id", "chr_name", "chrom_start"),
		filters    = c("snp_filter"),
		value      = snplist, 
		mart       = snpmart)
	names(snp_pos) <- c("snp", "chr", "position")
	snp_pos <- subset(snp_pos, select=c(snp, position))
	snp_pos <- subset(snp_pos, !duplicated(snp))
	sig$index <- 1:nrow(sig)
	sig <- merge(sig, snp_pos, by.x="snp1", by.y="snp", all.x=TRUE)
	sig <- merge(sig, snp_pos, by.x="snp2", by.y="snp", suff=c("1","2"), all.x=TRUE)
	sig <- sig[order(sig$index), ]
	sig <- subset(sig, select=-c(index))
	return(sig)
}


confInt <- function(n, alpha)
{
	k <- c(1:n)
	upper <- -log10(qbeta(alpha/2,k,n+1-k))
	lower <- -log10(qbeta((1-alpha/2),k,n+1-k))
	expect <- -log10((k-0.5)/n)
	return(data.frame(expect, lower, upper))
}


qqDat <- function(sig, alpha)
{
	a <- subset(sig, filter == 3)
	a <- a[order(a$pnest, decreasing=T), ]
	ci <- confInt(nrow(a), alpha)
	a <- data.frame(a, ci)

	b <- subset(sig, filter != 3)
	b <- b[order(b$pnest, decreasing=T), ]
	ci <- confInt(nrow(b), alpha)
	b <- data.frame(b, ci)
	return(rbind(b,a))
}


varianceComponents <- function(sig, geno, phen, loud=FALSE)
{
	a <- linearRegression(
		phen = phen[,sig$probeid],
		gen  = cbind(geno[,sig$pos1], geno[, sig$pos2]) + 1
	)
	if(loud) print(a)
	vars <- a$variances
	vars <- vars[names(vars) %in% c("aa", "ad", "da", "dd")]

	nom <- names(vars)[which.max(vars)]

	nom <- unlist(strsplit(toupper(nom), split=""))
	nom <- paste(nom, collapse=" x ")

	return(nom)
}


varianceComponentsBreakdown <- function(sig, geno, phen, loud=FALSE)
{
	a <- linearRegression(
		phen = phen[,sig$probeid],
		gen  = cbind(geno[,sig$pos1], geno[, sig$pos2]) + 1
	)
	if(loud) print(a)
	vars <- a$variances[-1]
	return(vars)
}


getVcBreakdown <- function(sig, geno, phen)
{
	l <- list()
	for(i in 1:nrow(sig))
	{
		cat(i, "of", nrow(sig), "\n")
		l[[i]] <- varianceComponentsBreakdown(sig[i,], xmat, resphen)
	}
	sig$vc <- l
	return(sig)
}


gpMaps <- function(sig, geno, phen)
{
	l <- list()
	m <- list()
	for(i in 1:nrow(sig))
	{
		cat(i, "of", nrow(sig), "\n")
		y <- phen[,sig$probeid[i]]
		x1 <- geno[,sig$pos1[i]]
		x2 <- geno[,sig$pos2[i]]
		l[[i]] <- tapply(y, list(x1, x2), function(x) mean(x, na.rm=T))
		m[[i]] <- table(x1, x2)
	}
	sig$gcm <- l
	sig$gcs <- m
	return(sig)
}


mergeBsgsRep <- function(bsgs, fehr, egcut)
{
	f <- subset(fehr, select=c(code, pfull, pnest, upper, vc, gcm, gcs))
	names(f) <- c("code", "pfull_fehr", "pnest_fehr", "upper_fehr", "vc_fehr", "gcm_fehr", "gcs_fehr")
	e <- subset(egcut, select=c(code, pfull, pnest, upper, vc, gcm, gcs))
	names(e) <- c("code", "pfull_egcut", "pnest_egcut", "upper_egcut", "vc_egcut", "gcm_egcut", "gcs_egcut")

	bsgs <- merge(bsgs, f, by="code", all.x=T)
	bsgs <- merge(bsgs, e, by="code", all.x=T)
	return(bsgs)
}


marginalSnpAssociations <- function(bsgs, marginal_list, threshold, probeinfo)
{
	m <- subset(marginal_list, pval < threshold)
	m <- m[order(m$pval, decreasing=FALSE), ]
	m <- subset(m, !duplicated(m$snp))

	i1 <- match(bsgs$snp1, m$snp)
	p1 <- m$probename[i1]
	g1 <- as.character(probeinfo$ILMN_GENE[match(p1, probeinfo$PROBE_ID)])
	bsgs$marginal_gene1 <- g1

	i2 <- match(bsgs$snp2, m$snp)
	p2 <- m$probename[i2]
	g2 <- as.character(probeinfo$ILMN_GENE[match(p2, probeinfo$PROBE_ID)])
	bsgs$marginal_gene2 <- g2

	return(bsgs)
}

# Load Data

load("~/repo/eQTL-2D/data/residuals_all.RData")
load("~/repo/eQTL-2D/data/clean_geno_final.RData")
load("~/repo/eQTL-2D/data/probeinfo_all.RData")
load("~/repo/eQTL-2D/filtering/marginal_lists/marginal_list.RData")
load("~/repo/eQTL-2D/replication/results/replication2_summarised.RData")
bsgs <- ReadOrig(
	"~/repo/eQTL-2D/replication/run/interactions_list2.RData",
	"sig",
	"BSGS"
)
egcut <- FormatRep(egcut)
egcut$set <- "EGCUT"
fehr <- FormatRep(fehr)
fehr$set <- "Fehrmann"

egcut <- qqDat(egcut, 0.05)
fehr <- qqDat(fehr, 0.05)
head(bsgs)
head(egcut)
head(fehr)

# Get VC
bsgs$probeid <- match(bsgs$probename, colnames(resphen))
bsgs <- getVcBreakdown(bsgs, xmat, resphen)
bsgs <- gpMaps(bsgs, xmat, resphen)

# Get associations for marginal SNPs
bsgs <- marginalSnpAssociations(bsgs, marginal_list, 1.29e-11, probeinfo_all)


# Get replication pvals
sig <- mergeBsgsRep(bsgs, fehr, egcut)
sig <- posData(sig)


# The VC lists are weird
for(i in 1:nrow(sig))
{
	sig$vc_fehr[[i]] <- sig$vc_fehr[[i]][[1]]
	sig$vc_egcut[[i]] <- sig$vc_egcut[[i]][[1]]
}


# Get sentinal SNPs and remove probes that tag the same genes
sig$code2 <- with(sig, paste(chr1, chr2, probegene))
sig <- sig[order(sig$code2, sig$pnest_egcut, decreasing=T), ]
sig <- subset(sig, !duplicated(code2) | filter == 3)


# Make sure that probechr, chr1 and chr2 are in 1:22
sig <- subset(sig, ((probechr %in% 1:22 & chr1 %in% 1:22 & chr2 %in% 1:22) | filter == 3))


# Make sure all SNPs are rsIDs
sig <- subset(sig, snp1 %in% grep("rs", snp1, value=TRUE) & snp2 %in% grep("rs", snp2, value=TRUE))
table(sig$filter != 3)


# Remove filter==3
sig_all <- sig
sig <- subset(sig, filter != 3)

sig_rep1 <- subset(sig, pnest_fehr > upper_fehr | pnest_egcut > upper_egcut)
sig_rep2 <- subset(sig, pnest_fehr > upper_fehr & pnest_egcut > upper_egcut)

save(sig, sig_all, sig_rep1, sig_rep2, file="~/repo/eQTL-2D/analysis/interaction_list_replication_summary.RData")
