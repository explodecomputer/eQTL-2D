# epi_empirical_analysis_fun.R

# analysis and interpretation of the output from epi_empirical.R
# joseph.powell@uq.edu.au



# provide summary formation for the output results
summarize.fun <- function(lf) {

	out <- matrix(0, nrow=length(lf), ncol=6)

	for(i in 1:length(lf)) {
		load(lf[i])
		index <- which(output$nclass==9 & output$minclass > 5 & output$LD < 0.1)
		foo <- output[index,]

		out[i,1] <- substr(lf[i], 1, 12)
		out[i,2] <- strsplit(lf[i], "_")[[1]][3]
		out[i,3] <- strsplit(lf[i], "_")[[1]][4]	

		out[i,4] <- nrow(foo)

		# Calculate lambda (median)
		Z <- qnorm(1-(foo$P/2))
		out[i,5] <- (median(na.omit(Z)))^2/0.456
		out[i,6] <- length(which(foo$P < 4.48e-6))

		# 

		rm(foo)
		rm(output)

		print(i)
	}

	out <- as.data.frame(out)
	names(out) <- c("probename", "snp1", "snp2", "nsnps", "lambda", "nthreshold")	
	return(out)
}
