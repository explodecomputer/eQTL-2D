#=======================================================#
#														#
#	epi_impute_analysis.R								#
#														#
#	analysis and plotting of the epi scan of the 		#
#	imputed blocks										#
#														#
#	joseph.powell@uq.edu.au		23/04/2013				#
#														#
#=======================================================#

# Impute output in: /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/epi_scan_imput_outputs
 
set3 <- load("/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/set3.RData")


#=======================================================#
#	pull in for information in the impute scan files	#
#=======================================================#

# Read in relavent imputed scan output (matched by probe and SNP 1 and SNP 2) 


impute_scan_filter.fun <- function(set3, dir) {
	# set3 : sig results from the genotyped SNPs
	# dir :	 dir of the impute scan outputs	


	out <- array(0, c(nrow(set3), 10))
	for (i in 1:nrow(set3)) {
		tmp <- read.table(paste(dir, "epi_impute_scan_", as.character(set3$probename[i]), "_", as.character(set3$snp1[i]), "_", as.character(set3$snp2[i]), ".txt", sep=""), header=T)

		# filter (rsq, ngenoclass, minclass size)
		index <- which(abs(tmp$rsq) < 0.1 & tmp$nclass==9 & tmp$minclass > 5)
		tmp <- tmp[index,]

		out[i,1:9] <- as.matrix(tmp[which.max(tmp$intP),])
		out[i,10] <- nrow(tmp)
		print(i)	

	}

	out <- as.data.frame(out)
	names(out) <- c(names(tmp), "ntests")	
	out <- cbind(set3, out)

	return(out)
}



impute_results <- impute_scan_filter.fun(set3, dir)
write.csv(umpute_results, "/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/impute_sig_results.csv", quote=F, row.names=F)








