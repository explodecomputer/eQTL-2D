ReadOrig <- function(filename, objname, setname)
{
	load(filename)
	sig <- get(objname)

	sig <- subset(sig, select=c(chr1, chr2, snp1, snp2, pos1, pos2, probename, probegene, probechr, pfull, pnest))
	sig$set <- setname
	sig$code <- with(sig, paste(probename, snp1, snp2))
	return(sig)
}


ReadRep <- function(filename, objname, setname)
{
	load(filename)
	sig <- get(objname)

	a <- sum(sig$replication_r^2 > 0.1)
	# cat(a, "pairs with high LD\n")

	b <- sum(sig$replication_nclass != 9)
	# cat(b, "pairs with < 9 classes\n")

	sig <- subset(sig, replication_r^2 < 0.1 & replication_nclass == 9, select=c(chr1, chr2, snp1, snp2, pos1, pos2, probename, probegene, probechr, replication_pfull, replication_pnest))

	names(sig) <- c("chr1", "chr2", "snp1", "snp2", "pos1", "pos2", "probename", "probegene", "probechr", "pfull", "pnest")
	sig$set <- setname
	sig$code <- with(sig, paste(probename, snp1, snp2))
	return(sig)
}


posData <- function(sig, bim)
{
	sig$position1 <- bim$V4[sig$pos1]
	sig$position2 <- bim$V4[sig$pos2]
	sig$snp1 <- as.character(sig$snp1)
	sig$snp2 <- as.character(sig$snp2)
	sig$probename <- as.character(sig$probename)
	sig$probegene <- as.character(sig$probegene)
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
	sig <- sig[order(sig$pnest, decreasing=T), ]
	ci <- confInt(nrow(sig), alpha)
	sig <- data.frame(sig, ci)
	return(sig)
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
	vars <- as.data.frame(as.list(a$variances))[,-1]
	sig <- cbind(sig, vars)

	return(sig)
}


getVcBreakdown <- function(sig, geno, phen)
{
	l <- list()
	for(i in 1:nrow(sig))
	{
		cat(i, "of", nrow(sig), "\n")
		l[[i]] <- varianceComponentsBreakdown(sig[i,], xmat, resphen)
	}
	l <- rbind.fill(l)
	return(l)
}


mergeBsgsRep <- function(bsgs, fehr, egcut)
{
	f <- subset(fehr, select=c(code, pfull, pnest, upper))
	names(f) <- c("code", "pfull_fehr", "pnest_fehr", "upper_fehr")
	e <- subset(egcut, select=c(code, pfull, pnest, upper))
	names(e) <- c("code", "pfull_egcut", "pnest_egcut", "upper_egcut")

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


load("~/repo/eQTL-2D/data/residuals_all.RData")
load("~/repo/eQTL-2D/data/clean_geno_final.RData")
bim <- read.table("~/repo/eQTL-2D/data/clean_geno_final.bim", colClasses=c("character", "character", "numeric", "numeric", "character", "character"))
load("~/repo/eQTL-2D/filtering/marginal_lists/marginal_list.RData")
bsgs <- ReadOrig(
	"~/repo/eQTL-2D/replication/run/interactions_list.RData",
	"sig",
	"BSGS"
)
egcut <- ReadRep(
	"~/repo/eQTL-2D/replication/results/EGCUT_replication.RData",
	"newsig",
	"EGCUT"
)
fehr <- ReadRep(
	"~/repo/eQTL-2D/replication/results/FehrmannHT12v3_replication.RData",
	"newsig",
	"Ferhmann"
)

bsgs <- posData(bsgs, bim)
fehr <- posData(fehr, bim)
egcut <- posData(egcut, bim)
egcut <- qqDat(egcut, 0.05)
fehr <- qqDat(fehr, 0.05)
head(bsgs)
head(egcut)
head(fehr)

# Get VC
bsgs$probeid <- match(bsgs$probename, colnames(resphen))
bsgs <- getVcBreakdown(bsgs, xmat, resphen)

# Get replication pvals
bsgs <- mergeBsgsRep(bsgs, fehr, egcut)

# Get associations for marginal SNPs
bsgs <- marginalSnpAssociations(bsgs, marginal_list, 1e-10, probeinfo)

sig <- bsgs

save(sig, file="~/repo/eQTL-2D/analysis/replication_summary.RData")




