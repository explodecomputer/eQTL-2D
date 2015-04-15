# epi_empirical_analysis_fun.R

# analysis and interpretation of the output from epi_empirical.R
# joseph.powell@uq.edu.au


##################################################################
##################################################################
##################################################################
# provide summary formation for the output results
summarize.fun <- function(lf) {

	out <- matrix(0, nrow=length(lf), ncol=11)

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
		out[i,5] <- round((median(na.omit(Z)))^2/0.456, 2)
		out[i,6] <- length(which(foo$P < 4.48e-6))
		# out[i,7] <- length(which(foo$P < 0.05/nrow(foo)))

		# Calculate N pairs with F_i > H0-F from 4.48x10-6
		F_thres <- qf(1-(4.48e-6), df1=4, df2=846)
		Q <- round(nrow(foo)*4.48e-6)

		if(Q==0) {

			F_sort <- sort(foo$F, decreasing=T)
			F_emp <- F_sort[1]
			P_emp <- round(-log10(1-pf(F_emp, df1=4, df2=846)),2)
			out[i,9] <- F_emp
			out[i,10] <- P_emp

			# Emp Type 1 error rate
			F_empN <- qf(1-(0.05/nrow(foo)), df1=4, df2=846)
			out[i,7] <- round(F_empN,2)
			out[i,8] <- length(which(F_sort > F_empN))

		}

		else{
			F_sort <- sort(foo$F, decreasing=T)
			F_emp <- F_sort[Q]
			P_emp <- round(-log10(1-pf(F_emp, df1=4, df2=842)),2)

			out[i,9] <- F_emp
			out[i,10] <- P_emp

			# Emp Type1 error
			F_empN <- qf(1-(0.05/nrow(foo)), df1=4, df2=846)
			out[i,7] <- round(F_empN,2)
			out[i,8] <- length(which(F_sort > F_empN))
		}
		
		# calculate the type 1 error rate
		F <- qf(1-(0.05), df1=4, df2=842) 
		out[i,11] <- round(length(which(foo$F > F))/nrow(foo),3)

			
		rm(foo)
		rm(output)

		print(i)
	}

	out <- as.data.frame(out)
	names(out) <- c("probename", "snp1", "snp2", "nsnps", "lambda", "nthreshold", "F_empNtests", "N_F_empNtests", "F_emp", "P_emp", "Type1")	
	return(out)
}




##################################################################
##################################################################
##################################################################
# Calculate the additive eQTL effects for each pair
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


##################################################################
##################################################################
##################################################################
# genome information summary
genome_sum.fun <- function(lf) {
	out <- matrix(0, nrow=length(lf), ncol=8)

	for(i in 1:length(lf)) {
		load(lf[i])

		out[i,1] <- substr(lf[i], 1, 12)
		out[i,2] <- strsplit(lf[i], "_")[[1]][3]
		out[i,3] <- strsplit(lf[i], "_")[[1]][4]	

		out[i,4] <- length(which(output$nclass == 9))
		out[i,5] <- length(which(output$minclass > 5))
		out[i,6] <- length(which(output$LD < 0.1))


		out[i,7] <- length(which(output$nclass==9 & output$minclass > 5 & output$LD < 0.1))
		out[i,8] <- nrow(output)

		rm(output)

		print(i)
	}

	out <- as.data.frame(out)
	names(out) <- c("probename", "snp1", "snp2", "nclass9", "minclass5", "LD01", "nsnpspass", "nsnps")	
	return(out)

}


##################################################################
##################################################################
##################################################################
# add information of the filter used and if it passed replication 
filter_add.fun <- function(sig, gs) {


	out <- matrix(0, nrow=nrow(gs), ncol=4)
	for(i in 1:nrow(gs)) {

		index <- which(sig$probename==gs$probename[i] & sig$snp1==gs$snp1[i] & sig$snp2==gs$snp2[i])
		if(length(index)!=1) {

			stop(print(i, "Error in the identification of matched probe and snps"))

		}

		out[i,1] <- sig$filter[index]
		out[i,2] <- sig$pnest_egcut[index]
		out[i,3] <- sig$pnest_fehr[index]
		out[i,4] <- sig$probegene[index]

	}

	out <- as.data.frame(out)
	names(out) <- c("filter", "pnest_egcut", "pnest_fehr", "gene")
	out <- cbind(gs, out)
	return(out)

}



##################################################################
##################################################################
##################################################################
# Calculate mean lambda for the multi-probe epi pairs

multi_lambda.fun <- function(gs, n) {
	# n = the number of epi pairs for a probe

	mp <- which(table(gs$probename)>n)

	out <- matrix(0, nrow=length(mp), ncol=6)
	for(i in 1:length(mp)) {

		foo <- gs[which(gs$probename==names(mp[i])),]
		out[i,5] <- round(mean(as.numeric(as.matrix(foo$lambda))),2)
		out[i,1:4] <- as.matrix(foo[1,c(1:3, 12)])
		out[i,6] <- nrow(foo)

	}

	out <- as.data.frame(out)
	names(out) <- c("probename", "snp1", "snp2", "gene", "meanlambda", "npairs")
	return(out)

}





##################################################################
##################################################################
##################################################################
# subset gs table to match the 30 paper significant pairs

type1_30.fun <- function(gs, sig30) {

	index <- rep(0, 30)
	gene <- 
	for(i in 1:nrow(sig30)) {
		index[i] <- which(gs$probename==as.character(sig30$Probe[i]) & gs$snp1==as.character(sig30$SNP1[i]) & gs$snp2==as.character(sig30$SNP2[i]))	
	}	

	gs30 <- gs[index,]
	gs30 <- cbind(sig30$GENE, gs30)
	names(gs30)[1] <- "gene"
	return(gs30)

}












