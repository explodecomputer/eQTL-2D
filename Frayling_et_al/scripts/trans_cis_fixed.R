# Trans mapping fixed for cis-snp
# Q: do we see an inflation in the p-vaues in a genome-wide analysis between 
# epistatsis snps when there are additive effects


# Run functions at the bottom of the script

# Read in data 
load("")			### Path dir to corrected expression data

# PLINK geno change function
ped <- read.table("")			# path to ped file 
map <- read.table("")			# path to map file

geno <- plink_to_012(ped, map)



# Models

trans_cis_epi.fun <- function(geno, snp1_id, probe, probe_id) {

	nsnps <- 				 ncol of geno data


	snp1 <- 			# fixed at the cis snp
	snp2 <- 			# this is to vary by i (1:nsnps)
	pheno <- 			# name / info for 1 of 2 probes			

	fullmod <- lm(pheno ~ as.factor(snp1) + as.factor(snp2) + as.factor(snp1):as.factor(snp2))
	redmod <- lm(pheno ~ as.factor(snp1) + as.factor(snp2))
	intmod <- anova(fullmod, redmod)





}



#=======================================================#
#		*******			FUNCTIONS 		********		#
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





