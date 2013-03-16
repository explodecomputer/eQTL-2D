#' Check files are present
#'
#' Makes sure that all the files required for the replication are present
#'
#' @param plink Path to \code{plink} executable
#' @param plinkfile Path to binary plinkfile (excluding any suffixes)
#' @param probefile Path to file with expression probe data
#' @param intlistfile Path to the \code{.RData} file that has all the target SNPs for replication
#'
#' @return Exits with an error if any files are missing
#' @export
CheckFiles <- function(plink, plinkfile, probefile, intlistfile)
{
	stopifnot(file.exists(plink))
	stopifnot(file.exists(plinkfile))
	stopifnot(file.exists(probefile))
	stopifnot(file.exists(intlistfile))
}


#' Read the probefile
#'
#' The probefile will be a plain text file, whitespace separated, with rows representing IDs and columns representing probes. The columns should be formatted as follows:
#' \describe{
#'  \item{\code{FID}}{Family ID (first row must be \code{FID})}
#'  \item{\code{IID}}{Individual ID column. All IID entries must be unique (first row must be \code{IID})}
#'  \item{\code{ILMN_xxxx}, \code{ILMN_xxxx} etc}{The rest of the columns are all the probes available in the set. Headers for the remaining columns must be the probe IDs.}}
#'
#' @param probefile Path to file with expression probe data
#'
#' @return Returns probe data as data.frame
#' @export
ReadProbeFile <- function(probefile)
{
	signal <- read.table(probefile, header=T)
	rownames(signal) <- signal$IID
	signal <- subset(signal, select=-c(IID, FID))
	return(signal)	
}


#' Load interaction list
#'
#' Loads the list of interactions (\code{.RData} file), and returns the subset which has SNPs and probes in common with the replication dataset
#'
#' @param intlistfile Path to the \code{.RData} file that has all the target SNPs for replication
#'
#' @return Returns \code{data.frame}
#' @export
LoadIntList <- function(intlistfile, plinkfile, probes)
{
	load(intlistfile)
	
	# Check overlap between SNPs in sig and replication data
	bim <- read.table(paste(plinkfile, ".bim", sep=""))
	repprobes <- colnames(probes)

	cat(sum(sig$probename %in% repprobes), "out of", nrow(sig), "interactions have expression probes in common with the replication data set\n")
	
	dim1 <- nrow(sig)
	sig1 <- subset(sig, snp1 %in% bim$V2 & snp2 %in% bim$V2)
	dim2 <- nrow(sig)

	cat(dim2, "out of", dim1, "interactions have SMPs in common with the replication data set\n")

	sig2 <- subset(sig1, probename %in% repprobes)
	return(sig2)
}


#' Read in SNPs required for replication
#'
#' Extracts all SNPs present in interaction list and reads in as 0/1/2 format matrix
#'
#' @param sig Output from \link{LoadIntList}
#' @param plink Path to \code{plink} executable
#' @param plinkfile Path to binary plinkfile (excluding any suffixes)
#'
#' @return Returns \code{matrix} of genotype data
#' @export
ExtractSNPs <- function(sig, plinkfile)
{
	require(snpStats)
	snplist <- with(sig, unique(c(as.character(snp1), as.character(snp2))))
	rawdata <- read.plink(bed=plinkfile, select.snps=snplist)
	geno <- apply(rawdata@.Data, 2, as.numeric)
	geno[geno==0] <- NA
	geno <- geno - 1

	fam <- read.table(paste(plinkfile, ".fam", sep=""), colClasses="character")
	bim <- read.table(paste(plinkfile, ".bim", sep=""), colClasses="character")
	rownames(geno) <- fam$V2
	colnames(geno) <- bim$V2

	return(geno)
}


#' Check probe and genotype data
#'
#' Checks that IIDs in geno and probes are matched, and reorders accordingly.
#'
#' @param probes Output from \link{ReadProbeFile}
#' @param geno Output from \link{ExtractSNPs}
#'
#' @return Returns \code{list} of geno \code{matrix} and probes \code{data.frame}
#' @export
DataChecks <- function(probes, geno)
{
	ids <- geno$ids
	snps <- geno$snps
	geno <- geno$geno

	# Check IDs are all unique
	stopifnot(!any(duplicated(out$ids)))

	# Check IDs are present in geno and probe data
	common_ids <- ids[ids %in% rownames(probes)]
	cat(length(common_ids), "individuals common in both geno and probe data\n")

	geno <- geno[rownames(geno) %in% common_ids, ]
	probes <- probes[rownames(probes) %in% common_ids, ]
	stopifnot(all(rownames(geno) %in% rownames(probes)))

	# reorder rows to match each other
	index <- match(rownames(probes), rownames(geno))
	geno <- geno[index, ]

	stopifnot(all(rownames(geno) == rownames(probes)))

	l <- list()
	l$geno <- geno
	l$ids <- rownames(geno)
	l$snps <- colnames(geno)
	l$probes <- probes

	return(l)
}


#' Replication statistical tests
#'
#' Calculates allele frequencies, correlation between SNPs, class sizes
#'
#' @param geno \code{matrix} of genotype data
#' @param probes \code{data.frame} of probes
#' @param sig Output from \link{LoadIntList}
#' @param i Which row of \code{sig} to run the analysis on
#'
#' @return Returns \code{data.frame} with row \code{i} complete
#' @export
ReplicationTests <- function(geno, probes, sig, i)
{
	# Extract data
	snp1 <- geno[colnames(geno) == sig$snp1]
	snp2 <- geno[colnames(geno) == sig$snp2]
	probe <- probes[colnames(probes) == sig$probename]

	# Summary statistics
	tab <- table(snp1 + 3*snp2)

	sig$replication_p1[i] <- mean(snp1) / 2
	sig$replication_p2[i] <- mean(snp2) / 2
	sig$replication_r[i] <- cor(snp1, snp2, use="pair")
	sig$replication_nclass[i] <- length(tab)
	sig$replication_minclass[i] <- min(tab)
	sig$replication_nid[i] <- sum(!is.na(snp1) & !is.na(snp2))

	# Statistical tests
	fullmod <- lm(probe ~ as.factor(snp1) * as.factor(snp2))
	margmod <- lm(probe ~ as.factor(snp1) + as.factor(snp2))
	fulltest <- summary(fullmod)$fstatistic
	inttest <- anova(margmod, fullmod)

	sig$replication_pfull[i] <- -log10(pf(fulltest[1], fulltest[2], fulltest[3], low=FALSE))
	sig$replication_pnest[i] <- -log10(inttest$P[2])

	return(sig)
}


#' Run replication analysis
#'
#' Test all SNP pairs in interaction list in replication dataset.
#'
#' @param sig Output from \link{LoadIntList}
#' @param checked Output from \link{DataChecks}
#'
#' @return Returns \code{data.frame} with new columns for results from replication data
#' @export
RunReplication <- function(sig, checked)
{
	sig$replication_pfull <- NA
	sig$replication_pnest <- NA
	sig$replication_p1 <- NA
	sig$replication_p2 <- NA
	sig$replication_r <- NA
	sig$replication_nclass <- NA
	sig$replication_minclass <- NA
	sig$replication_nid <- NA

	geno <- checked$geno
	probes <- checked$probes

	for(i in 1:nrow(sig))
	{
		cat(i, "of", nrow(sig), "\n")
		sig <- ReplicationTests(snp1, snp2, probe, sig, i)
	}
	return(sig)
}
