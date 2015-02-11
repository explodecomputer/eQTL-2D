# epi_empirical_analysis.R

# analysis and interpretation of the output from epi_empirical.R
# joseph.powell@uq.edu.au

# set the data dir
setwd("/Users/jpowell/repo/eQTL-2D/investigation/data/output/")




lf <- list.files()

summarize.fun <- function(lf) {

	out <- matrix(0, nrow=length(lf), ncol=4)

	for(i in 1:length(lf)) {
		load(lf[i])
		index <- which(output$nclass==9 & output$minclass > 5 & output$LD < 0.1)
		foo <- output[index,]

		# Calculate lambda (median)

		out[i,1] <- substr(lf[i], 1, 12)
		out[i,2] <- nrow(foo)
		Z <- qnorm(1-(foo$P/2))
		out[i,3] <- (median(na.omit(Z)))^2/0.456
		out[i,4] <- length(which(foo$P < 4.48e-6))

		rm(foo)
		rm(output)

		print(i)
	}

	out <- as.data.frame(out)
	names(out) <- c("probename", "nsnps", "lambda", "nthreshold")	
	return(out)
}


out <- summarize.fun(lf)




# Calculate the additive eQTL effects for each pair


load("/Users/jpowell/repo/eQTL-2D/investigation/data/")

add_cal.fun <- function(sig, bsgs, geno, bim) {

	out <- matrix(0, nrow=nrow(sig), ncol=6)
	for(i in 1:nrow(sig)) {

		probe <- as.character(sig$probename[i])
		snp1 <- as.character(sig$snp1[i])
		snp2 <- as.character(sig$snp2[i])

		pheno <- bsgs[,which(names(bsgs)==probe)]
		geno1 <- geno[,which(names(geno)==snp1)]
		geno2 <- geno[,which(names(geno)==snp2)]

		# check everything is the correct length
		if(length(geno1)!=846 | length(geno2)!=846 | length(pheno)!=846) {
			stop(print(paste(probe, " incorrect data length")))
		}

		# Run single marker additive model
		out[i,4] <- summary(lm(pheno~geno1))$coefficients[2,4]
		out[i,5] <- summary(lm(pheno~geno2))$coefficients[2,4]
		out[i,6] <- length(which(out[i,4:5] < 0.0000001))
		# add names
		out[i,1] <- probe
		out[i,2] <- snp1	
		out[i,3] <- snp2	

	}

	out <- as.data.frame(out)
	names(out) <- c("probe", "snp1", "snp2", "addpval1", "addpval2", "nadd")
	return(out)


}

add_test <- add_cal.fun(sig, bsgs, geno, bim)

