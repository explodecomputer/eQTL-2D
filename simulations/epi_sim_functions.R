#=======================================================#
#-------------------------------------------------------#
#														#
#	epi_sim_functions.R									#
#														#
#	Set of functions to run a simulation on imputerd 	#
#	regions. 											#
#														#
#	joseph.powell@uq.edu.au								#
#	V1. March.2013										#
#														#
#-------------------------------------------------------#
#=======================================================#


#=======================================================#
#	CHOOSE THE PAIRWISE COMBINATION OF SNPS TO TEST		#
#=======================================================#

snp_pick.fun  <- function(
	type,			# type of epi imteraction; cis-cis, cis-trans, trans-trans
	probe_info,		# dataframe of the probe info
	geno_info) 		# genotype snp information
	{


	# select a probe
	probe <- probe_info[sample(c(1:nrow(probe_info)), 1),]

		
	# select SNPs according to type
	if(type=="cis_cis") {
		# make the geno snp selection panel
		index <- which(geno_info$chr==probe$CHROMOSOME_NEW & geno_info$position > probe$PROBE_START-1000000 & geno_info$position < probe$PROBE_START+1000000)
		snp_selection <- geno_info[index,]

		snps_out <- snp_selection[sample(c(1:nrow(snp_selection)), 2, replace=F),]
		out <- list(snps_out, probe)
		return(out)
	}

	if(type=="cis_trans") {
		# make the geno snp selection panel
		index <- which(geno_info$chr==probe$CHROMOSOME_NEW & geno_info$position > probe$PROBE_START-1000000 & geno_info$position < probe$PROBE_START+1000000)
		snp_selection_cis <- geno_info[index,]
		snp_selection_trans <- geno_info[-index,]

		snps_out_cis <- snp_selection_cis[sample(c(1:nrow(snp_selection_cis)), 1),]
		snps_out_trans <- snp_selection_trans[sample(c(1:nrow(snp_selection_trans)), 1),]

		snps_out <- rbind(snps_out_cis, snps_out_trans)
		out <- list(snps_out, probe)
		return(out)
		
	}

	if(type=="trans_trans") {
		# make the geno snp selection panel
		index <- which(geno_info$chr==probe$CHROMOSOME_NEW & geno_info$position > probe$PROBE_START-1000000 & geno_info$position < probe$PROBE_START+1000000)
		snp_selection <- geno_info[-index,]

		snps_out <- snp_selection[sample(c(1:nrow(snp_selection)), 2, replace=F),]
		out <- list(snps_out, probe)
		return(out)

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
#	.FUN FOR THE SNP BY SNP PAIRWISE 4DF AND 8DF MODELS #
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
		print(i)
	}

	out <- as.data.frame(out)
	names(out) <- c("snp1", "snp2", "rsq", "df1", "df2", "fullP", "intP", "nclass", "minclass")
	return(out)

}



#=======================================================#
#	TEST FOR THE ASSOCIATION BETWEN THE GENOTYPED SNSP	#
#=======================================================#



epi_geno.fun <- function(
	geno_snp1, 		# genotypes for SNP1
	geno_snp2,		# genotypes for SNP2
	probe, 			# matched probe phenotype 
	snp1,			# names 
	snp2 			# names
	) {

	out <- array("NA", c(1, 9))

	# Information
	out[1,1] <- snp1
	out[1,2] <- snp2

	# rsq
	out[1,3] <- round(cor(geno_snp1, geno_snp2), 4)

	# 4 and 8df tests	

	fullmod <- lm(probe ~ as.factor(geno_snp1) + as.factor(geno_snp2) + as.factor(geno_snp1):as.factor(geno_snp2))
	redmod <- lm(probe ~ as.factor(geno_snp1) + as.factor(geno_snp2))
	# This is the interaction terms on their own (nested test)
	intmod <- anova(redmod, fullmod)	

	# Extract statistics	
	tmp <- summary(fullmod)$fstatistic
	out[1,4] <- tmp[2]
	out[1,5] <- tmp[3]
	out[1,6] <- round(-log10(pf(tmp[1], tmp[2], tmp[3], low=F)),4)		

	out[1,7] <- round(-log10(intmod$Pr[2]), 4)

	# class sizes
	out[1,8] <- length(table(geno_snp1 + 3*geno_snp2))
	out[1,9] <- min(table(geno_snp1 + 3*geno_snp2))

	out <- as.data.frame(out)
	names(out) <- c("snp1", "snp2", "rsq", "df1", "df2", "fullP", "intP", "nclass", "minclass")
	return(out)

}








