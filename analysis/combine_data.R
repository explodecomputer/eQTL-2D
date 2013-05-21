# Make list of SNPs with
SNP position
SNP chromosome
SNP gene
original interaction tests
replication interaction tests
previous associations
probe gene


replicationOverlap <- function(bsgs, fehr, egcut)
{
	# Make a variable that indicates level of replication
	# - not present
	# - no replication
	# - 5% FDR in one
	# - 5% FDR in two
	# - Bonferroni in one
	# - Bonferroni in two

	egcut$bonf <- FALSE
	egcut$bonf[egcut$pnest > -log10(0.05/nrow(egcut))] <- TRUE
	fehr$bonf <- FALSE
	fehr$bonf[fehr$pnest > -log10(0.05/nrow(fehr))] <- TRUE

	egcut$fdr <- FALSE
	egcut$fdr[egcut$pnest > egcut$upper] <- TRUE
	fehr$fdr <- FALSE
	fehr$fdr[fehr$pnest > fehr$upper] <- TRUE


	bsgs$rep <- "None"
	bsgs$sets <- 0

	# FDR in one

	bsgs$rep[bsgs$code %in% egcut$code[egcut$fdr] | bsgs$code %in% fehr$code[fehr$fdr]] <- "FDR"
	bsgs$sets[bsgs$code %in% egcut$code[egcut$fdr] | bsgs$code %in% fehr$code[fehr$fdr]] <- 1
	bsgs$sets[bsgs$code %in% egcut$code[egcut$fdr] & bsgs$code %in% fehr$code[fehr$fdr]] <- 2

	bsgs$rep[bsgs$code %in% egcut$code[egcut$bonf] | bsgs$code %in% fehr$code[fehr$bonf]] <- "Bonf"
	bsgs$sets[bsgs$code %in% egcut$code[egcut$bonf] | bsgs$code %in% fehr$code[fehr$bonf]] <- 1
	bsgs$sets[bsgs$code %in% egcut$code[egcut$bonf] & bsgs$code %in% fehr$code[fehr$bonf]] <- 2

	table(bsgs$rep, bsgs$sets)
	return(bsgs)
}

newAssoc <- function(sig, marginal_list)
{
	sig$code1 <- paste(sig$snp1, sig$probename)
	sig$code2 <- paste(sig$snp2, sig$probename)
	marginal_list$code <- paste(marginal_list$snp, marginal_list$probename)
	sig$known1 <- "known"
	sig$known1[! sig$code1 %in% marginal_list$code] <- "new"
	sig$known2 <- "known"
	sig$known2[! sig$code2 %in% marginal_list$code] <- "new"

	sig <- subset(sig, select=-c(code1, code2))
	return(sig)
}



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


mergeBsgeRep <- function(bsgs, egcut, fehr)
{

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

# Get 





bsgs <- replicationOverlap(bsgs, fehr, egcut)
sig <- probesWithRep(bsgs)



sig <- newAssoc(sig, marginal_list)
sig <- subset(sig, select=c(snp1, chr1, position1, pos1, known1, snp2, chr2, position2, pos2, known2, probename, probegene, probechr, pfull, pnest, rep, sets, a., d., .a, .d, aa, ad, da, dd))

save(sig, file="~/repo/eQTL-2D/analysis/interaction_list_summary.RData")




