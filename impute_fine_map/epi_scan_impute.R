#=======================================================#
#-------------------------------------------------------#
#														#
#	epi_scan_impute.R									#
#														#
#	Scan for 4 and 8df models (epistasis) on the 	 	#
#	Imputed blocks										#
#														#
#	joseph.powell@uq.edu.au								#
#	V1. Feb.2013										#
#														#
#-------------------------------------------------------#
#=======================================================#



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
#				READ IN THE DATA AND ARGS 				#
#=======================================================#

n <- as.numeric(commandArgs(T)[1])
pheno <-  commandArgs(T)[2]
set3 <- commandArgs(T)[3]
blockdir <- commandArgs(T)[4]
outdir <- commandArgs(T)[5]

load(pheno)
load(set3)


#=======================================================#
#		GET THE PHENOTYPE AND GENOTYPE DATA FOR N		#
#=======================================================#

ped1 <- read.table(paste(blockdir, as.character(set3$snp1[n]), ".ped", sep=""), header=F)
ped2 <- read.table(paste(blockdir, as.character(set3$snp2[n]), ".ped", sep=""), header=F)
 
map1 <- read.table(paste(blockdir, as.character(set3$snp1[n]), ".map", sep=""), header=F)
map2 <- read.table(paste(blockdir, as.character(set3$snp2[n]), ".map", sep=""), header=F)



block1 <- plink_to_012.fun(ped1, map1)
block2 <- plink_to_012.fun(ped2, map2)

pheno <- resphen[,which(colnames(resphen)==as.character(set3$probename[n]))]


#=======================================================#
#					RUN EPI_SCAN.FUN 			 		#
#=======================================================#

epi_scan_out <- epi_scan.fun(block1, block2, pheno)


#=======================================================#
#		WRITE OUT THE OUTPUT INTOTHE OUTDIR 	 		#
#=======================================================#

write.table(epi_scan_out, paste(outdir, "epi_impute_scan_", as.character(set3$probename[n]), "_",  as.character(set3$snp1[n]), "_", as.character(set3$snp2[n]), ".txt", sep=""), quote=F, row.names=F)


