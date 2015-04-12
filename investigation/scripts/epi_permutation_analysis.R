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

	setwd(dir)
	lf <- list.files(dir)
	out <- array(0, c(length(lf), 6))

	for(i in 1:length(lf)) {
		load(lf[i])

		probe <- substr(lf[i], 1, 12)
		snp1 <- strsplit(lf[i], "_")[[1]][[3]]
		snp2 <- strsplit(lf[i], "_")[[1]][[4]]

		# Filter the output	
		index <- which(output$nclass == 9 & output$minclass > 4 & output$LD < 0.01)
		foo <- output[index, ]
		foo <- foo[order(foo$P),]
		foo$P <- -log10(foo$P)

		pair <- which(sig$probename==probe & sig$snp1==snp1 & sig$snp2==snp2)
		pemp <- which.min(abs(foo$P-sig$pnest[pair]))/10000000

		out[i,1] <- sig$probename[pair]
		out[i,2] <- sig$probegene[pair]
		out[i,3] <- sig$snp1[pair]
		out[i,4] <- sig$snp2[pair]
		out[i,5] <- round(sig$pnest[pair], 3)
		out[i,6] <- round(-log10(pemp), 3)

		print(i)
	}

	out <- as.data.frame(out)
	names(out) <- c("probename", "probegene", "snp1", "snp2", "pnest", "pemp")
	return(out)
}


pemp_out <- analysis.fun("~/repo/eQTL-2D/investigation/data/output_permutation/", sig)
write.table(pemp_out, "~/repo/eQTL-2D/investigation/data/pemp_out.txt", quote=F, row.names=F)


png(filename="~/repo/eQTL-2D/investigation/docs/figures/pemp.png", width=600, height=600)
hist(as.numeric(as.matrix(pemp_out$pemp)), breaks=25, 
	xlab="permutation p - value",
	main="501 Pairs 10,000,000")
dev.off()

