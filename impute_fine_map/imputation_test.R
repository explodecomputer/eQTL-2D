#=======================================================#
#-------------------------------------------------------#
#														#
#	imputation_test.R									#
#														#
#	extract genotype regions for the significant pairs	#
#	convert the SNP formats and test the full regions	#
#														#
#	joseph.powell@uq.edu.au								#
#	V1. Feb.2013										#
#														#
#-------------------------------------------------------#
#=======================================================#

#=======================================================#
#		READ IN THE FILTERED DATASETS (FROM GIB)		#
#=======================================================#


load("/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/agg_filtered.RData")
snp_info <- read.table("/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Var_eQTL/bsgs_imputed_R2_80_cleaned_stage2_chr_all_SNP_info.txt", header=T)
snp_info_geno <- read.table("/ibscratch/wrayvisscher/josephP/BSGS/Genotype/GWAS.map", header=T)
load("/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/residuals.RData")



#=======================================================#
#	IDENTIFY A SINGLE SNP PAIR FOR EACH INTERACTION		#
#=======================================================#


SNP_pair_find.fun <- function(
	data 		# significant hits
	) {

	probes <- as.character(unique(data$probename))
	out <- array(0, c(length(probes), ncol(data)))

	for(i in 1:length(probes)) {
		tmp <- 	data[which(data$probename==probes[i]), ]

		out[i,] <- as.matrix(tmp[which.max(tmp$pnest),])
	}

	out <- as.data.frame(out)
	names(out) <- names(data)
	return(out)	
}

#===# Use

#snp_pair <- SNP_pair_find.fun(set2)


#=======================================================#
#	EXTRACT THE RELAVENT SNP REGIONS FROM IMPUTED		#
#=======================================================#

make_snp_region_files.fun <- function(
	snp, 		# SNP pairs files
	) {

	for(i in 1:nrow(snp)) {

		# Data surrounding SNP 1
		snp1 <- as.character(snp$snp1[i])	
		system(paste("/clusterdata/apps/plink/plink-1.07-x86_64/plink --bfile /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Clean_R2_Imputed_BSGS_data/bsgs_imputed_R2_80_cleaned_stage2_chr_all --snp ", snp1, " --window 100 --recode12 --out /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/snp_region_data/", snp1,  sep=""))	

		# Data surrounding SNP 2
		snp2 <- as.character(snp$snp2[i])	
		system(paste("/clusterdata/apps/plink/plink-1.07-x86_64/plink --bfile /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Clean_R2_Imputed_BSGS_data/bsgs_imputed_R2_80_cleaned_stage2_chr_all --snp ", snp2, " --window 100 --recode12 --out /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/snp_region_data/", snp2,  sep=""))	

		# Merge the two files
		p_name <- as.character(snp$probename[i])
		system(paste("/clusterdata/apps/plink/plink-1.07-x86_64/plink --file /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/snp_region_data/", snp1, " --merge /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/snp_region_data/", snp2, ".ped /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/snp_region_data/", snp2, ".map", " --recode --out /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/snp_region_data/", snp1, "_", snp2, "_", p_name, sep=""))	

	}
}



#=======================================================#
#	IDENTIFYALL SNPS NOT IN THE IMPUTE DATA AND PROXYS	#
#=======================================================#



missing_snps.fun <- function(
	lf, 		# List of files in the output directory	
	set2,		# the set of epistasis associations (filtered)

	info1,		# snp infor from genotyped (only) data
	info2) {		# snp info from the imputed data

	tmp1 <- NULL

	for(i in 1:nrow(set2)) {
		i1 <- which(lf==paste(as.character(set2$snp1[i]), ".ped", sep=""))
		if(length(i1)==0) {
			tmp1 <- rbind(tmp1, as.character(set2$snp1[i]))	
		}
	}

	for(i in 1:nrow(set2)) {
		i2 <- which(lf==paste(as.character(set2$snp2[i]), ".ped", sep=""))
		if(length(i2)==0) {
			tmp1 <- rbind(tmp1, as.character(set2$snp2[i]))	
		}
	}


	# find proxy
	snps <- unique(tmp1)
	out <- array(NA, c(length(snps), 2))

	for(k in 1:length(snps)) {

		index <- snp_info_geno[which(snp_info_geno$rs_id==snps[k]),]

		if(index$chr !< 23)
	}

	# mainly on the X chromosome - function unfinished

}

#=======================================================#
#	IDENTIFY ALL IMPUTED (FULL) EPI ASSOCIATIONS LEFT	#
#=======================================================#



epi_impute_match.fun <- function(
	lf, 		# List of files in the output directory	
	set2		# the set of epistasis associations (filtered)
	) {		# snp info from the imputed data

	out <- array("NA", c(nrow(set2)))

	for(i in 1:nrow(set2)) {

		i1 <- which(lf==paste(as.character(set2$snp1[i]), ".ped", sep=""))
		i2 <- which(lf==paste(as.character(set2$snp2[i]), ".ped", sep=""))

		if(length(i1)!=0 & length(i2)!=0) {
			out[i] <- "GOOD"
		}

		else{
			out[i] <- "BAD"
		}

	}


	set3 <- cbind(set2, out)
	set3 <- as.data.frame(set3)
	names(set3)[ncol(set3)] <- "IMPUTE_OK"
	return(set3)
}


lf <- list.files("/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/snp_region_data/")
set3 <- epi_impute_match.fun(lf, set2)
save(set3, file="set3.RData")


#=======================================================#
#		CONVERT THE PED FORMT TO A 0, 1, 2, FORMAT 		#
#=======================================================#


plink_to_012.fun <- function(
	pedfile, 		# Address of ped file
	mapfile) {		# Address of map file

	# Read files in
	pformat <- read.table(pedfile, header=FALSE, colClasses="character")
	map <- read.table(mapfile, header=FALSE, colClasses=c("character", "character", "numeric", "numeric"))

	ids <- pformat[, 1:6]
	nid <- nrow(ids)
	pformat <- pformat[, -c(1:6)]
	index <- seq(1, ncol(pformat), 2)
	geno <- matrix(0, nid, length(index))

	# Convert to 0, 1, 2 format
	for(i in 1:length(index)) {
		snp <- pformat[,c(index[i], index[i]+1)]
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








