# epi_empirical_analysis_fun.R

# analysis and interpretation of the output from epi_empirical.R
# joseph.powell@uq.edu.au



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
# lambda GC check 

gc_check.fun <- function(lf, sig, bsgs, geno, bim) {
	
	out <- matrix(0, nrow=length(lf), ncol=8)

	for(i in 1:length(lf)) {
	

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


		fullmod <- lm(pheno ~ as.factor(geno1) + as.factor(geno2) + as.factor(geno1):as.factor(geno2))
		redmod <- lm(pheno ~ as.factor(geno1) + as.factor(geno2))
		intmod <- anova(fullmod, redmod)

		load(lf[i])
		index <- which(output$nclass==9 & output$minclass > 5 & output$LD < 0.1)
		foo <- output[index,]

		out[i,1] <- substr(lf[i], 1, 12)
		out[i,2] <- strsplit(lf[i], "_")[[1]][3]
		out[i,3] <- strsplit(lf[i], "_")[[1]][4]	

		out[i,4] <- nrow(foo)

		# Calculate lambda (median)
		Z1 <- qnorm(1-(foo$P/2))
		F4 <- qf(1-foo$P, df1=4, df2=846)

		lambdaC <- round((median(na.omit(Z1)))^2/0.456, 2)
		lambdaF <- round((median(na.omit(F4)))/0.84, 2)
		out[i,5] <- lambdaC
		out[i,6] <- lambdaF

		F_bonf <- intmod$F[2]
		out[i,7] <- 1-pchisq(qchisq(pf(F_bonf, df1=4, df2=842), 1)/lambdaC, 1)	

		# calculate the P from the adjusted F lambda		
		out[i,8] <- 1-pf(F_bonf/lambdaF, df1=4, df2=842)			

		print(i)
	}
		
	out <- as.data.frame(out)
	names(out) <- c("probename", "snp1", "snp2", "nsnps", "lambdaC", "lambdaF", "PlamC", "PlamF")	
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


	out <- matrix(0, nrow=nrow(gs), ncol=9)
	for(i in 1:nrow(gs)) {

		index <- which(sig$probename==gs$probename[i] & sig$snp1==gs$snp1[i] & sig$snp2==gs$snp2[i])
		if(length(index)!=1) {

			stop(print(i, "Error in the identification of matched probe and snps"))

		}

		out[i,1] <- sig$filter[index]
		out[i,2] <- sig$pnest_egcut[index]
		out[i,3] <- sig$pnest_fehr[index]
		out[i,4] <- sig$probegene[index]
		out[i,5] <- sig$chr1[index]
		out[i,6] <- sig$chr2[index]
		out[i,7] <- sig$pos1[index]
		out[i,8] <- sig$pos2[index]
		out[i,9] <- sig$probechr[index]

	}

	out <- as.data.frame(out)
	names(out) <- c("filter", "pnest_egcut", "pnest_fehr", "gene", "chr1", "chr2", "pos1", "pos2", "probechr")
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
	for(i in 1:nrow(sig30)) {
		index[i] <- which(gs$probename==as.character(sig30$Probe[i]) & gs$snp1==as.character(sig30$SNP1[i]) & gs$snp2==as.character(sig30$SNP2[i]))	
	}	

	gs30 <- gs[index,]
	gs30 <- cbind(sig30$GENE, gs30)
	names(gs30)[1] <- "gene"
	return(gs30)

}



#=======================================================#
#	RUN THE ANALYSIS FOR STEP ONE OF THE INVESTIGATION 	#
#=======================================================#


analysis.fun <- function(probe, snp1, snp2, bsgs, geno, bim) {


	pheno <- bsgs[,which(names(bsgs)==probe)]
	geno1 <- geno[,which(names(geno)==snp1)]
	geno2 <- geno[,which(names(geno)==snp2)]

	# check everything is the correct length
	if(length(geno1)!=846 | length(geno2)!=846 | length(pheno)!=846) {
		stop(print(paste(i, " incorrect data length")))
	}

	# Run single marker additive model
	add1 <- summary(lm(pheno~geno1))$coefficients[2,4]
	add2 <- summary(lm(pheno~geno2))$coefficients[2,4]

	# perform genome-wide anaysis
	p <- which.min(c(add1, add2))
	if(p==1) {
		snp_fix <- snp1
		snp_other <- snp2
	}
	if(p==2) {
		snp_fix <- snp2
		snp_other <- snp1
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


	out <- matrix(0, nrow=ncol(g), ncol=4)
	telliter <- 1000
	for(k in 1:nrow(out)) {

		fit <- summary(lm(pheno~g[,k]))
		out[k,1] <- fit$coefficients[2,4]
		out[k,2] <- fit$df[2]
		out[k,3] <- fit$coefficients[2,2]
		out[k,4] <- fit$coefficients[2,3]

		if(k %% telliter==0){
			print(k)
		}

	}

	out <- as.data.frame(out)
	names(out) <- c("pval", "df", "se", "fstat")
	return(out)

}



#=======================================================#
#		CONVERT THE PED FORMT TO A 0, 1, 2, FORMAT 		#
#=======================================================#

plink_to_012.fun <- function(
	ped, 		# ped file
	map) {		# map file

	ids <- ped[, 1:6]
	nid <- nrow(ids)
	ped <- ped[, -c(1:6)]
	index <- seq(1, ncol(ped), 2)
	geno <- matrix(0, nid, length(index))

	# Convert to 0, 1, 2 format
	for(i in 1:length(index)) {
		snp <- ped[,c(index[i], index[i]+1)]
		x <- array(NA, nid)
		snp[snp == "0"] <- NA

		i0 <- snp[,1] == 1 & snp[,2] == 1
		i2 <- snp[,1] == 2 & snp[,2] == 2
		i1 <- (snp[,1] == 1 & snp[,2] == 2) | (snp[,1] == 2 & snp[,2] == 1)
		x[i0] <- 0
		x[i1] <- 1
		x[i2] <- 2
		geno[, i] <- x
	}

	colnames(geno) <- map$V2
	rownames(geno) <- ids$V2
	return(geno)
}




TableOfTruth.fun <- function(data, sig30) {


	datasig <- data[which(as.numeric(as.matrix(data$PlamF)) < 4.48e-6 & data$pemp < 4.48e-6),]

	# convert the p's to non -log 10 scale
	datasig$pnest_egcut <- 10^-as.numeric(as.matrix(datasig$pnest_egcut))
	datasig$pnest_fehr <- 10^-as.numeric(as.matrix(datasig$pnest_fehr))
	
	# Calculate F's
	datasig$egcutF <- qf(1-datasig$pnest_egcut, df1=4, df2=891)
	datasig$fehrF <- qf(1-datasig$pnest_fehr, df1=4, df2=1240)
	
	# Calculate F adj for lambdaFgc
	datasig$egcutFgc <- datasig$egcutF/as.numeric(as.matrix(datasig$lambdaF))
	datasig$fehrFgc <- datasig$fehrF/as.numeric(as.matrix(datasig$lambdaF))
	
	# Calculate adjusted p's
	datasig$pFgc_egcut <- 1-pf(datasig$egcutFgc, df1=4, df2=891)
	datasig$pFgc_fehr <- 1-pf(datasig$fehrFgc, df1=4, df2=1240)

	# set any 0's to 10e-20
	datasig$pFgc_fehr[which(datasig$pFgc_fehr==0)] <- 1.0e-25
	datasig$pFgc_egcut[which(datasig$pFgc_egcut==0)] <- 1.0e-25


	# Calculate the fisher meta-pvalue
    fit <- fisher.comb.fun(datasig[,37:38])
	datasig$combPFgc <- round(-log10(1-pchisq(fit$S, 2)), 1)

	# add cis-trans etc
	datasig$cis_trans <- "cis_trans"
	cc <- which(datasig$chr1==datasig$chr2)
	datasig$cis_trans[cc] <- "cis_cis"
	tt <- which(datasig$chr1!=datasig$probechr & datasig$chr2!=datasig$probechr)
	datasig$cis_trans[tt] <- "trans_trans"

	# Add Y/N for table 1 of hemani
	datasig$table1 <- "No"
	for(i in 1:nrow(sig30)) {
		index <- which(as.character(datasig$probename)==as.character(sig30$Probe[i]) & as.character(datasig$snp1)==as.character(sig30$SNP1[i]) & as.character(datasig$snp2)==as.character(sig30$SNP2[i]))
		if(length(index)==1) {

			datasig$table1[index] <- "Yes"
		}
	}	

	# sort by replication p-value
	datasig <- datasig[ order(datasig$combPFgc, decreasing=T),]
	return(datasig)	

}



#####################
#
# Fisher 1948 combined test of significance for independent tests
# Mosteller F, Fisher RA. Combining independent tests of significance.
# The American Statistician, Vol. 2, No. 5 (Oct., 1948), pp. 30-31.
#
# 'aka' the standard fisher combined method
#
#####################
fisher.comb.fun <- function (pvals, method = c("fisher"), p.corr = c("bonferroni", 
    "BH", "none"), zero.sub = 1e-05, na.rm = FALSE, mc.cores = NULL) 
{
    stopifnot(method %in% c("fisher"))
    stopifnot(p.corr %in% c("none", "bonferroni", "BH"))
    stopifnot(all(pvals >= 0, na.rm = TRUE) & all(pvals <= 1, 
        na.rm = TRUE))
    stopifnot(zero.sub >= 0 & zero.sub <= 1 || length(zero.sub) != 
        1)
    if (is.null(dim(pvals))) 
        stop("pvals must have a dim attribute")
    p.corr <- ifelse(length(p.corr) != 1, "BH", p.corr)
    pvals[pvals == 0] <- zero.sub
    if (is.null(mc.cores)) {
        fisher.sums <- data.frame(do.call(rbind, apply(pvals, 
            1, f.sum.fun, zero.sub = zero.sub, na.rm = na.rm)))
    }
    else {
        fisher.sums <- mclapply(1:nrow(pvals), function(i) {
            f.sum.fun(pvals[i, ], zero.sub = zero.sub, na.rm = na.rm)
        }, mc.cores = mc.cores)
        fisher.sums <- data.frame(do.call(rbind, fisher.sums))
    }
    rownames(fisher.sums) <- rownames(pvals)
    fisher.sums$p.value <- 1 - pchisq(fisher.sums$S, df = 2 * 
        fisher.sums$num.p)
    fisher.sums$p.adj <- switch(p.corr, bonferroni = p.adjust(fisher.sums$p.value, 
        "bonferroni"), BH = p.adjust(fisher.sums$p.value, "BH"), 
        none = fisher.sums$p.value)
    return(fisher.sums)
}


# sum the fisher test statistic
f.sum.fun <- function (p, zero.sub = 1e-05, na.rm = FALSE) 
{
    if (any(p > 1, na.rm = TRUE) || any(p < 0, na.rm = TRUE)) 
        stop("You provided bad p-values")
    stopifnot(zero.sub >= 0 & zero.sub <= 1 || length(zero.sub) != 
        1)
    p[p == 0] <- zero.sub
    if (na.rm) 
    p <- p[!is.na(p)]

    S = -2 * sum(log(p))
    res <- data.frame(S = S, num.p = length(p))
    return(res)
}


# Function to perform the standard fisher analysis on k pvalues, looped over n 'genes'
fisher_combine_analysis.fun <- function(pval) {

    out <- array(0, c(nrow(pval), 2))

    for(i in 1:nrow(out)) {

        fit <- fisher.comb.fun(as.matrix(pval[i,]))
        fit2 <- f.sum.fun(fit$p.adj)    
        out[i,1] <- round(fit2$S, 3)
        out[i,2] <- round(-log10(1-pchisq(fit2$S, ncol(pval))), 3)
    }


    out <- as.data.frame(out)
    names(out) <- c("Combined_stat", "Combined_-log10pval")
    return(out)
}




