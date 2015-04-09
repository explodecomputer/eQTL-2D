# epi_permutation_analysis.R
# code and functions to perform the analysis of the permutation output for each of the stage 2 probes
# joseph.powell@uq.edu.au
# V1: 09/04/2015

# Read in the data for the 501 probe / snp pairs
load('~/repo/eQTL-2D/investigation/data/sig_501_data.RData')





# Analysis function
analysis.fun <- function(dir, sig) {
	# dir: the directory where the *RData files are located
	# sig: data on the significant permutation results

	lf <- list.files(dir)
	out <- array(0, c(length(lf), ))



	for(i in 1:length(lf)) {





		print(i)
	}


}

