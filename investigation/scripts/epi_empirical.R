# author: joseph.powell@uq.edu.au
# epi_empirical.R

# analysis to determine the empirical threshold

# Read in data 
# sig data
load('~/repo/eQTL-2D/investigation/data/sig_501_data.RData')	

# Genotype data
load('~/repo/eQTL-2D/data/geno.RData')
geno <- as.data.frame(geno)

# Probe (residuals) data
load('~/repo/eQTL-2D/data/bsgs_egcut_fehr_data.RData')
bsgs <- as.data.frame(phenlist[[1]])
rm(phenlist)






for(i in 1:nrow(sig)) {

	probe <- as.character(sig$probename[i])
	snp1 <- as.character(sig$snp1[i])
	snp2 <- as.character(sig$snp2[i])

	probe <- bsgs[,which(names(bsgs)==probe)]
	geno1 <- geno[,which(names(geno)==snp1)]
	geno2 <- geno[,which(names(geno)==snp2)]

	# check everything is the correct length
	if(length(geno1)!=846 | length(geno2)!=846 | length(probe)!=846) {
		stop(print(paste(i, " incorrect data length")))
	}

	# Run single marker additive model
	add1 <- summary(lm(probe~geno1))$coefficients[2,4]
	add2 <- summary(lm(probe~geno2))$coefficients[2,4]

	# perform genome-wide anaysis
	p <- which.min(c(add1, add2))
	if(p==1) {
		snp_fix <- snp1
		snp_other <- snp2
	}
	if(p==2) {
		snp_fix <- snp2
		snp_other <- snp1
	}

	# make the 'remaining' geno
	i1 <- which(bim$V1==bim$V1[which(bim$V2==snp_other)])
	MB_plus5 <- bim$V4[which(bim$V2==snp_fix)]+5000000 
	MB_minus5 <- bim$V4[which(bim$V2==snp_fix)]-5000000 
	chr <- bim$V1[which(bim$V2==snp_fix)]
	i2 <- which(bim$V1==chr & bim$V4 > MB_minus5 & bim$V4 < MB_plus5)
	index <- unique(c(i1,i2))
	
	# none conflicted geno matrix
	g <- geno[,-index]


	out <- matrix(0, nrow=ncol(g), ncol=6)
	for(k in 1:nrow(out)) {

		fit <- summary(lm(probe~g[,k]))
		out[k,1] <- fit$coefficients[2,4]

	}


}
