# epi_permutation.R
# code and functions to perform the permutation analysis for each of the stage 1 probes
# joseph.powell@uq.edu.au
# V1: 22/02/2015



args        		<- commandArgs(T)
data   				<- args[1]			
iter				<- as.numeric(args[2])
outdir     			<- args[3]	

load(data)

probe <- as.character(sig$probename[iter])
snp1 <- as.character(sig$snp1[iter])
snp2 <- as.character(sig$snp2[iter])


perm_analysis.fun <- function(probe, snp1, snp2, bsgs, geno, bim) {

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



	out <- matrix(0, nrow=10000000, ncol=5)
	telliter <- 10000
	for(k in 1:nrow(out)) {

		g <- sample(geno_other)
		fullmod <- lm(pheno ~ as.factor(geno_fix) + as.factor(g) + as.factor(geno_fix):as.factor(g))
		redmod <- lm(pheno ~ as.factor(geno_fix) + as.factor(g))
		intmod <- anova(fullmod, redmod)

		out[k,1] <- round(intmod$F[2],2)		# store F-statistics
		out[k,2] <- intmod$Pr[2]				# store P-values

		out[k,3] <- length(table(geno_fix + 3*g))   # nclass size
		out[k,4] <- min(table(geno_fix + 3*g))		# min class size

		out[k,5] <- round(cor(geno_fix, g, use="pairwise.complete.obs"),2)^2		# LD / correlation between 2 snps 

		if(k %% telliter==0){
			print(k)
		}

	}

	out <- as.data.frame(out)
	names(out) <- c("F", "P", "nclass", "minclass", "LD")
	return(out)

}



output <- analysis.fun(probe, snp1, snp2, bsgs, geno, bim)
save(output, file=paste(outdir, probe, "_", snp1, "_", snp2, "_permutation_output.RData", sep=""))






