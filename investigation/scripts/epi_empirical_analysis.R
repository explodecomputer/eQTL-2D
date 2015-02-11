# epi_empirical_analysis.R

# analysis and interpretation of the output from epi_empirical.R
# joseph.powell@uq.edu.au

# set the data dir
setwd("/Users/jpowell/repo/eQTL-2D/investigation/data/output/")




lf <- list.files()

summarize.fun <- function(lf) {

	out <- matrix(0, nrow=length(lf), ncol=6)

	for(i in 1:length(lf)) {
		load(lf[i])
		index <- which(output$nclass==9 & output$minclass > 5 & output$LD < 0.1)
		foo <- output[index,]

		# Calculate lambda (median)

		out[i,1] <- substr(lf[i], 1, 12)
		out[i,2] <- nrow(foo)
		out[i,3] <- min(foo$P) 	
		

		rm(foo)
		rm(output)
	}


Z <- qnorm(1-(foo$P/2))
lambda <- (median(na.omit(Z)))^2/0.456

}

for