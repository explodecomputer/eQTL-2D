# Trans mapping fixed for cis-snp
# Q: do we see an inflation in the p-vaues in a genome-wide analysis between 
# epistatsis snps when there are additive effects


# Run functions at the bottom of the script

# Read in data 
load("/Users/jpowell/repo/eQTL-2D/data/bsgs_egcut_fehr_data.RData")			### Path dir to corrected expression data
pheno <- phenlist[[1]]

# PLINK geno change function
load("/Users/jpowell/repo/eQTL-2D/data/geno.RData")




#=======================================================#
#		RUN ANALYSIS USING FUNCTIONS LISTED BELOW		#
#=======================================================#

# TMEM149 - ILMN_1786426
# CIS SNPS - rs8106959 

TMEM149_out <- trans_cis_epi.fun(geno, "rs8106959", pheno, "ILMN_1786426")
MBLN1_out <- trans_cis_epi.fun(geno, "rs13069559", pheno, "ILMN_2313158")

write.csv(TMEM149_out, "/Users/jpowell/repo/eQTL-2D/Frayling_et_al/data_files/TMEM149_out.csv", quote=F, row.names=F)
write.csv(MBLN1_out, "/Users/jpowell/repo/eQTL-2D/Frayling_et_al/data_files/MBLN1_out.csv", quote=F, row.names=F)




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



#=======================================================#
#		ANALYSIS OF EPI SCAN BY FIXED CIS SNP 			#
#=======================================================#

trans_cis_epi.fun <- function(geno, snp1_id, probe, probe_id) {


	snp1 <- geno[,which(colnames(geno)==snp1_id)]
	g <- geno[,-which(colnames(geno)==snp1_id)]				
	probe <- pheno[,which(colnames(pheno)==probe_id)]						

	nsnps <- ncol(g)	# ncol of geno data

	out <- array(0, c(nsnps, 5))				# blank array of results
	telliter <- 1000

	for(i in 1:nsnps) {
	
		snp2 <- g[,i]			# this is to vary by i (1:nsnps)
	

		fullmod <- lm(probe ~ as.factor(snp1) + as.factor(snp2) + as.factor(snp1):as.factor(snp2))
		redmod <- lm(probe ~ as.factor(snp1) + as.factor(snp2))
		intmod <- anova(fullmod, redmod)

		out[i,1] <- round(intmod$F[2],2)		# store F-statistics
		out[i,2] <- -log10(intmod$Pr[2])		# store P-values

		out[i,3] <- length(table(snp1 + 3*snp2))   	# nclass size
		out[i,4] <- min(table(snp1 + 3*snp2))		# min class size

		out[i,5] <- round(cor(snp1, snp2, use="pairwise.complete.obs"),2)^2		# LD / correlation between 2 snps 


		# print iteraction
		if(i %% telliter==0) {
			cat(paste("iteraction ", i, " complete\n"))
		}


	}

	out <- as.data.frame(out)
	names(out) <- c("F", "P", "nclass", "minclass", "LD")
	return(out)

}



