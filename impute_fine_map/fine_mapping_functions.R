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
	info, 		# SNP pairs files
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
		print(i)
	}

	out <- as.data.frame(out)
	names(out) <- c("snp1", "snp2", "rsq", "df1", "df2", "fullP", "intP", "nclass", "minclass")
	return(out)

}





#=======================================================#
#	COMPARE THE IMPUTED RESULTS AGAINST GENOTYPED		#
#=======================================================#

impute_vs_geno.fun <- function(
	set3, n) {		#	set3 data
	

	out <- NULL

	for(i in 1:nrow(set3)) {

		epi <- read.table(paste("/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/epi_scan_imput_outputs/epi_impute_scan_", as.character(set3$probename[i]), "_", as.character(set3$snp1[i]), "_", as.character(set3$snp2[i]), ".txt", sep=""), header=T)
		index <- which(epi$snp1==as.character(set3$snp1[i]) & epi$snp2==as.character(set3$snp2[i]))


		if(length(index)!=1) {

			print("Hello, what's going on with this pair?")

		}

		else {

			tmp <- epi[index,]	
			out <- rbind(out, tmp)

		}
			# Filter
#			index <- which(epi$nclass==9 & epi$minclass>4)
#			epi <- epi[index,]
#			out <- rbind(out, epi[which.max(epi$intP),])
		print(i)
	
		}


		return(out)
}





impute_vs_geno.fun <- function(
	set3, n) {		#	set3 data
	

	out <- NULL
	out2 <- NULL

	for(i in 1:nrow(set3)) {

		epi <- read.table(paste("/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/epi_scan_imput_outputs/epi_impute_scan_", as.character(set3$probename[i]), "_", as.character(set3$snp1[i]), "_", as.character(set3$snp2[i]), ".txt", sep=""), header=T)
		# Filter
		index <- which(epi$nclass==9 & epi$minclass>4)
		epi <- epi[index,]
		out <- rbind(out, epi[which.max(epi$intP),])
		out2 <- c(out2, nrow(epi))

		print(i)
	
		}

		return(out)


}








