#=======================================================#
#-------------------------------------------------------#
#														#
#	fine_mapping_functions.R							#
#														#
#	run the analysis to fine map the interactions for 	#
#	the list of probes and snps given in the wood et al	#
#	list.												#
#														#
#	joseph.powell@uq.edu.au								#
#	V1. March.2014										#
#														#
#-------------------------------------------------------#
#=======================================================#

#=======================================================#
#	EXTRACT THE RELAVENT SNP REGIONS FROM IMPUTED		#
#=======================================================#

make_snp_region_files.fun <- function(
	info 		# SNP pairs files
	) {

	for(i in 3:nrow(info)) {

		# Data surrounding SNP 1
		snp1 <- as.character(info$SNP1[i])	
		system(paste("/clusterdata/apps/plink/plink-1.07-x86_64/plink --bfile /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Clean_R2_Imputed_BSGS_data/bsgs_imputed_R2_80_cleaned_stage2_chr_all --snp ", snp1, " --window 10 --recode12 --out /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/snp_region_data/", snp1,  sep=""))	

		# Data surrounding SNP 2
		snp2 <- as.character(info$SNP2[i])	
		system(paste("/clusterdata/apps/plink/plink-1.07-x86_64/plink --bfile /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Clean_R2_Imputed_BSGS_data/bsgs_imputed_R2_80_cleaned_stage2_chr_all --snp ", snp2, " --window 10 --recode12 --out /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/snp_region_data/", snp2,  sep=""))	

		# Merge the two files
		#p_name <- as.character(info$Probe[i])
		#system(paste("/clusterdata/apps/plink/plink-1.07-x86_64/plink --file /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/snp_region_data/", snp1, " --merge /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/snp_region_data/", snp2, ".ped /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/snp_region_data/", snp2, ".map", " --recode --out /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/snp_region_data/", snp1, "_", snp2, "_", p_name, sep=""))	

	}
}




#=======================================================#
#		CONVERT THE PED FORMT TO A 0, 1, 2, FORMAT 		#
#=======================================================#

plink_to_012.fun <- function(
	ped, 		# ped file
	map) {		# map file

	ids <- ped[, 1:6]
	nid <- nrow(ids)
	ped <- ped[, -c(1:6)]
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




#=======================================================#
#	RUN THE SNP BY SNP PAIRWISE 4DF AND 8DF MODELS 		#
#=======================================================#


epi_scan.fun <- function(
	block1, 		# SNP1 block
	block2,			# SNP2 block
	probe 			# matched probe phenotype 
	) {

	# Check the sample ids match
	out <- array(NA, c(ncol(block1)*ncol(block2), 9))
	c <- 0

	for(i in 1:ncol(block1)) {
		snpi <- block1[,i]

		for(k in 1:ncol(block2)) {
			c <- c+1
			snpk <- block2[,k]	

			# check the SNP names are different
			if(colnames(block1)[i]==colnames(block2)[k]) {

			#	print("matching snp ids")	
				out[c,] <- "NA" 	
			}

			else {
				# Information
				out[c,1] <- colnames(block1)[i]
				out[c,2] <- colnames(block2)[k]

				# rsq
				out[c,3] <- round(cor(snpi, snpk), 4)

				# 4 and 8df tests	

				fullmod <- lm(probe ~ as.factor(snpi) + as.factor(snpk) + as.factor(snpi):as.factor(snpk))
#				out[i,10] <- summary(fullmod)$r.squared
				redmod <- lm(probe ~ as.factor(snpi) + as.factor(snpk))
				# This is the interaction terms on their own (nested test)
				intmod <- anova(redmod, fullmod)	

				# Extract statistics	
				tmp <- summary(fullmod)$fstatistic
				out[c,4] <- tmp[2]
				out[c,5] <- tmp[3]
				out[c,6] <- round(-log10(pf(tmp[1], tmp[2], tmp[3], low=F)),4)		

				out[c,7] <- round(-log10(intmod$Pr[2]), 4)

				# class sizes
				out[c,8] <- length(table(snpi + 3*snpk))
				out[c,9] <- min(table(snpi + 3*snpk))


			}	

		}	
		#print(i)
	}

	out <- as.data.frame(out)
	names(out) <- c("snp1", "snp2", "rsq", "df1", "df2", "fullP", "intP", "nclass", "minclass")
	return(out)

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


