# author: joseph.powell@uq.edu.au
# epi_empirical_regress.R

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


analysis_regress.fun <- function(probe, snp1, snp2, bsgs, geno, bim) {

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
		geno_other <- geno2
	}
	if(p==2) {
		snp_fix <- snp2
		snp_other <- snp1
		geno_fix <- geno2
		geno_other <- geno1
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


	# make adjusted phenotype (by 1 snp)
	pheno_adj1 <- array(NA, length(pheno))
	index1 <- !is.na(geno_fix)
	fit1 <- lm(pheno~geno_fix)
	pheno_adj1[index1] <- fit1$residuals	

	# make adjusted phenotype (by snps)
	pheno_adj2 <- array(NA, length(pheno))
	fit2 <- lm(pheno~geno_fix+geno_other)
	pheno_adj2[index1] <- fit2$residuals


	out <- matrix(0, nrow=ncol(g), ncol=7)
	telliter <- 100
	for(k in 1:50000){#nrow(out)) {

		# run test for adjusted phenotype 1
		snp2 <- g[,k]			# this is to vary by k (1:nsnps)
		fullmod1 <- lm(pheno_adj1 ~ as.factor(geno_fix) + as.factor(snp2) + as.factor(geno_fix):as.factor(snp2))
		redmod1 <- lm(pheno_adj1 ~ as.factor(geno_fix) + as.factor(snp2))
		intmod1 <- anova(fullmod1, redmod1)

		# run test for adjusted phenotype 2
		fullmod2 <- lm(pheno_adj2 ~ as.factor(geno_fix) + as.factor(snp2) + as.factor(geno_fix):as.factor(snp2))
		redmod2 <- lm(pheno_adj2 ~ as.factor(geno_fix) + as.factor(snp2))
		intmod2 <- anova(fullmod2, redmod2)

		out[k,1] <- round(intmod1$F[2],2)		# store F-statistics
		out[k,2] <- intmod1$Pr[2]				# store P-values

		out[k,3] <- round(intmod2$F[2],2)		# store F-statistics
		out[k,4] <- intmod2$Pr[2]				# store P-values


		out[k,5] <- length(table(geno_fix + 3*snp2))   	# nclass size
		out[k,6] <- min(table(geno_fix + 3*snp2))		# min class size

		out[k,7] <- round(cor(geno_fix, snp2, use="pairwise.complete.obs"),2)^2		# LD / correlation between 2 snps 

		if(k %% telliter==0){
			print(k)
		}

	}

	out <- as.data.frame(out)
	names(out) <- c("F1", "P1", "F2", "P2", "nclass", "minclass", "LD")
	return(out)

}


output <- analysis_regress.fun(probe, snp1, snp2, bsgs, geno, bim)
save(output, file=paste(outdir, probe, "_", snp1, "_", snp2, "_regress_output.RData", sep=""))







