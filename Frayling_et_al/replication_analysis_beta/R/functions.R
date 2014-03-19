# 2. Alter the subsequent functions to analyse the new data format
# 3. add the functions to run the following;
# 3a. Predicting genotypes
# 3b. Fitting the effects before and after the inc snp
# 3c. Where possible the analysis of the inc snp and the 'other' snp
# 3d. anythink else?

#' Check files are present
#'
#' Makes sure that all the files required for the replication are present
#'
#' @param plinkfile Path to binary plinkfile (excluding any suffixes)
#' @param probefile Path to file with expression probe data
#' @param intlistfile Path to the \code{.RData} file that has all the target SNPs for replication
#'
#' @return Exits with an error if any files are missing
#' @export
CheckFiles <- function(plinkfile, probefile, intlistfile)
{
	stopifnot(file.exists(paste(plinkfile, ".map", sep="")))
	stopifnot(file.exists(paste(plinkfile, ".ped", sep="")))
	stopifnot(file.exists(probefile))
	stopifnot(file.exists(intlistfile))
}

#' Read in and convert genotype file
#'
#' @param plinkfile Path to ped and map file (excluding any suffixes)
#'
#' @return Exits with an error if any files are missing
#' @export
GenoIN <- function(plinkfile) 
{
	ped <- read.table(paste(plinkfile, ".ped", sep=""), header=F)
	map <- read.table(paste(plinkfile, ".map", sep=""), header=F)

	nsnps <- nrow(map)
	n <- ncol(ped)-(2*nsnps)

	ids <- ped[, 1:n]
	nid <- nrow(ids)
	ped <- ped[, -c(1:n)]
	index <- seq(1, ncol(ped), 2)
	geno <- matrix(0, nid, length(index))

	# Convert to 0, 1, 2 format
	for(i in 1:length(index)) {
		snp <- ped[,c(index[i], index[i]+1)]
		x <- array(NA, nid)
		snp[snp == "0"] <- NA

		i0 <- snp[,1] == 1 & snp[,2] == 1
		i2 <- snp[,1] == 2 & snp[,2] == 2
		i1 <- (snp[,1] == 1 & snp[,2] == 2) | (snp[,1] == 2 & snp[,2] == 1)
		x[i0] <- 0
		x[i1] <- 1
		x[i2] <- 2
		geno[, i] <- x
	}

	colnames(geno) <- map$V2
	rownames(geno) <- ids$V2
	return(geno)
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
	
	repprobes <- colnames(probes)

	cat(sum(sig$probename %in% repprobes), "out of", nrow(sig), "interactions have expression probes in common with the replication data set\n")
	
	dim1 <- nrow(sig)
	sig1 <- subset(sig, snp1 %in% colnames(geno) & snp2 %in% colnames(geno))
	dim2 <- nrow(sig1)

	cat(dim2, "out of", dim1, "interactions have SNPs in common with the replication data set\n")

	sig2 <- subset(sig1, probename %in% repprobes)
	return(sig2)
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
	ids <- rownames(geno)
	snps <- colnames(geno)

	# Check IDs are all unique
	stopifnot(!any(duplicated(ids)))

	# Check IDs are present in geno and probe data
	common_ids <- ids[ids %in% rownames(probes)]
	cat(length(common_ids), "individuals common in both geno and probe data\n")

	geno <- geno[rownames(geno) %in% common_ids, ]
	probes <- probes[rownames(probes) %in% common_ids, ]
	stopifnot(all(rownames(geno) %in% rownames(probes)))

	# reorder rows to match each other
	index <- match(rownames(probes), rownames(geno))
	geno <- geno[index, ]
	ids <- ids[index]

	stopifnot(all(rownames(geno) == rownames(probes)))
	stopifnot(all(rownames(geno) == ids))

	l <- list()
	l$geno <- geno
	l$ids <- rownames(geno)
	l$snps <- colnames(geno)
	l$probes <- probes

	return(l)
}


#' Replication statistical tests
#'
#' Calculates allele frequencies, correlation between SNPs, class sizes, variance components
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
	require(noia)
	# Extract data
	snp1 <- geno[, colnames(geno) == sig$snp1[i]]
	snp2 <- geno[, colnames(geno) == sig$snp2[i]]
	probe <- probes[, colnames(probes) == sig$probename[i]]

	# Summary statistics
	tab <- table(snp1 + 3*snp2)
	gcm <- tapply(probe, list(snp1, snp2), function(x) { mean(x, na.rm=T)})
	gcs <- table(snp1, snp2)
	mod <- linearRegression(probe, cbind(snp1, snp2)+1)

	sig$replication_p1[i] <- mean(snp1, na.rm=T) / 2
	sig$replication_p2[i] <- mean(snp2, na.rm=T) / 2
	sig$replication_r[i] <- cor(snp1, snp2, use="pair")
	sig$replication_nclass[i] <- length(tab)
	sig$replication_minclass[i] <- min(tab, na.rm=T)
	sig$replication_nid[i] <- sum(!is.na(snp1) & !is.na(snp2))

	# Statistical tests
	fullmod <- lm(probe ~ as.factor(snp1) * as.factor(snp2))
	margmod <- lm(probe ~ as.factor(snp1) + as.factor(snp2))
	fulltest <- summary(fullmod)$fstatistic
	inttest <- anova(margmod, fullmod)

	sig$replication_pfull[i] <- -log10(pf(fulltest[1], fulltest[2], fulltest[3], low=FALSE))
	sig$replication_pnest[i] <- -log10(inttest$P[2])

	l <- list()
	l$sig <- sig
	l$gcm <- gcm
	l$gcs <- gcs
	l$mod <- mod

	return(l)
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

	gcm <- list()
	gcs <- list()
	mod <- list()

	for(i in 1:nrow(sig))
	{
		cat(i, "of", nrow(sig), "\n")
		out <- ReplicationTests(geno, probes, sig, i)
		sig <- out$sig
		gcm[[i]] <- out$gcm
		gcs[[i]] <- out$gcs
		mod[[i]] <- out$mod
	}

	l <- list()
	l$sig <- sig
	l$gcm <- gcm
	l$gcs <- gcs
	l$mod <- mod

	return(l)
}



#' Run genotype prediction
#'
#' predict the genogtypes of one snp from the other two.
#'
#' @param sig Output from \link{LoadIntList}
#' @param checked Output from \link{DataChecks}
#'
#' @return Returns \code{data.frame} with new columns for results from replication data
#' @export
GenoPrediction <- function(sig, checked)
{
	
	sig_new <- sig[which(sig$pfull==0),]

	out <- array(0, c(nrow(sig_new)/2, 4))
	k <- 0
	for(i in seq(1, nrow(sig_new), 2)) {
		k <- k+1
		y <- which(checked$snps==as.character(sig_new$snp2[i]))
		x1 <- which(checked$snps==as.character(sig_new$snp1[i]))
		x2 <- which(checked$snps==as.character(sig_new$snp1[i+1]))

		if(length(y)==1 & length(x1)==1 & length(x2)==1) {

			out[k,1] <- summary(lm(checked$geno[,y] ~ as.factor(checked$geno[,x1])))$adj.r.squared
			out[k,2] <- summary(lm(checked$geno[,y] ~ as.factor(checked$geno[,x2])))$adj.r.squared
			out[k,3] <- summary(lm(checked$geno[,y] ~ as.factor(checked$geno[,x1]) + as.factor(checked$geno[,x2])))$adj.r.squared
			out[k,4] <- summary(lm(checked$geno[,y] ~ as.factor(checked$geno[,x1]) + as.factor(checked$geno[,x2]) + as.factor(checked$geno[,x1]):as.factor(checked$geno[,x2])))$adj.r.squared

		}
	}
	out <- round(out, 2)	
	return(out)
}


#' Run before and after analysis
#'
#' estimate the effects before and after correcting for the inc seq snp
#'
#' @param sig Output from \link{LoadIntList}
#' @param checked Output from \link{DataChecks}
#'
#' @return Returns \code{data.frame} with new columns for results from replication data
#' @export
CorrectionTest <- function(sig, checked)
{
	require(noia)
	sig_new <- sig[which(sig$pfull==0),]

	geno <- checked$geno
	probes <- checked$probes

	gcm <- list()
	gcs <- list()
	mod <- list()

	k <- 0
	out2 <- list()

	for(i in seq(1, nrow(sig_new), 2)) {
		k <- k+1
		y <- which(checked$snps==as.character(sig_new$snp2[i]))
		x1 <- which(checked$snps==as.character(sig_new$snp1[i]))
		x2 <- which(checked$snps==as.character(sig_new$snp1[i+1]))
		p1 <- which(names(checked$probes)==as.character(sig_new$probename[i]))

		if(length(y)==1 & length(x1)==1 & length(x2)==1 & length(p1)==1) {

			out <- NULL

			out$snp1 <- as.character(sig_new$snp1[i])
			out$snp2 <- as.character(sig_new$snp1[i+1])
			out$incsnp <- as.character(sig_new$snp2[i])
			out$probe <- as.character(sig_new$probename[i])

			padj <- summary(lm(probes[,p1] ~ as.factor(geno[,y])))$residuals
			porj <- probes[,p1] 

			# Extract data
			snp1 <- geno[, x1]
			snp2 <- geno[, x2]

			# Summary statistics
			tab <- table(snp1 + 3*snp2)
			out$tab <- tab
			out$gcm1 <- tapply(padj, list(snp1, snp2), function(x) { mean(x, na.rm=T)})
			out$gcm2 <- tapply(porj, list(snp1, snp2), function(x) { mean(x, na.rm=T)})

			out$gcs <- table(snp1, snp2)
			out$mod1 <- linearRegression(padj, cbind(snp1, snp2)+1)
			out$mod2 <- linearRegression(porj, cbind(snp1, snp2)+1)
			
			# find a way to store the output
			out$replication_p1 <- mean(snp1, na.rm=T) / 2
			out$replication_p2 <- mean(snp2, na.rm=T) / 2
			out$replication_r <- cor(snp1, snp2, use="pair")
			out$replication_nclass <- length(tab)
			out$replication_minclass <- min(tab, na.rm=T)
			out$replication_nid <- sum(!is.na(snp1) & !is.na(snp2))

			# Statistical tests 1
			fullmod1 <- lm(padj ~ as.factor(snp1) * as.factor(snp2))
			margmod1 <- lm(padj ~ as.factor(snp1) + as.factor(snp2))
			fulltest1 <- summary(fullmod1)$fstatistic
			inttest1 <- anova(margmod1, fullmod1)

			# Statistical tests 2
			fullmod2 <- lm(porj ~ as.factor(snp1) * as.factor(snp2))
			margmod2 <- lm(porj ~ as.factor(snp1) + as.factor(snp2))
			fulltest2 <- summary(fullmod2)$fstatistic
			inttest2 <- anova(margmod2, fullmod2)

			# find a way to store the output

			out$replication_pfull <- -log10(pf(fulltest1[1], fulltest1[2], fulltest1[3], low=FALSE))
			out$replication_pnest <- -log10(inttest1$P[2])

			out$replication_pfull2 <- -log10(pf(fulltest2[1], fulltest2[2], fulltest2[3], low=FALSE))
			out$replication_pnest2 <- -log10(inttest2$P[2])
		}

	out2[[k]] <- out 		

	}

	return(out2)

}


