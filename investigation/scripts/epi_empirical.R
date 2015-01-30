# author: joseph.powell@uq.edu.au
# epi_empirical.R

# analysis to determine the empirical threshold

# Read in data 
# sig data
#load('~/repo/eQTL-2D/investigation/data/sig_501_data.RData')	

# Genotype data
#load('~/repo/eQTL-2D/data/geno.RData')
#geno <- as.data.frame(geno)

# Probe (residuals) data
#load('~/repo/eQTL-2D/data/bsgs_egcut_fehr_data.RData')
#bsgs <- as.data.frame(phenlist[[1]])
#rm(phenlist)
# save(bim, bsgs, fam, freq, geno, sig, file="~/repo/eQTL-2D/investigation/data/investigation_data.RData")

# To run, e.g.:
# cd analysis/run
# R --no-save --args /path/to/investigation/RData array_iteration /path/to/output < epi_empirical.R

args        		<- commandArgs(T)
data   				<- args[1]			
iter				<- as.numeric(args[2])
outdir     			<- args[3]	

load(data)

probe <- as.character(sig$probename[iter])
snp1 <- as.character(sig$snp1[iter])
snp2 <- as.character(sig$snp2[iter])


analysis.fun <- function(probe, snp1, snp2, bsgs, geno, bim) {

	pheno <- bsgs[,which(names(bsgs)==probe)]
	geno1 <- geno[,which(names(geno)==snp1)]
	geno2 <- geno[,which(names(geno)==snp2)]

	# check everything is the correct length
	if(length(geno1)!=846 | length(geno2)!=846 | length(pheno)!=846) {
		stop(print(paste(probe, " incorrect data length")))
	}

	# Run single marker additive model
	add1 <- summary(lm(pheno~geno1))$coefficients[2,4]
	add2 <- summary(lm(pheno~geno2))$coefficients[2,4]

	# perform genome-wide anaysis
	p <- which.min(c(add1, add2))
	if(p==1) {
		snp_fix <- snp1
		snp_other <- snp2
		geno_fix <- geno1
	}
	if(p==2) {
		snp_fix <- snp2
		snp_other <- snp1
		geno_fix <- geno2
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


	out <- matrix(0, nrow=ncol(g), ncol=5)
	telliter <- 100
	for(k in 1:nrow(out)) {

		snp2 <- g[,k]			# this is to vary by k (1:nsnps)
		fullmod <- lm(pheno ~ as.factor(geno_fix) + as.factor(snp2) + as.factor(geno_fix):as.factor(snp2))
		redmod <- lm(pheno ~ as.factor(geno_fix) + as.factor(snp2))
		intmod <- anova(fullmod, redmod)

		out[k,1] <- round(intmod$F[2],2)		# store F-statistics
		out[k,2] <- intmod$Pr[2]				# store P-values

		out[k,3] <- length(table(geno_fix + 3*snp2))   	# nclass size
		out[k,4] <- min(table(geno_fix + 3*snp2))		# min class size

		out[k,5] <- round(cor(geno_fix, snp2, use="pairwise.complete.obs"),2)^2		# LD / correlation between 2 snps 

		if(k %% telliter==0){
			print(k)
		}

	}

	out <- as.data.frame(out)
	names(out) <- c("F", "P", "nclass", "minclass", "LD")
	return(out)

}


output <- analysis.fun(probe, snp1, snp2, bsgs, geno, bim)
save(output, file=paste(outdir, probe, "_output.RData", sep=""))







