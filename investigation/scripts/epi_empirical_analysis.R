# epi_empirical_analysis.R

# analysis and interpretation of the output from epi_empirical.R
# joseph.powell@uq.edu.au

library(xtable)

# read in the functions
source("/Users/jpowell/repo/eQTL-2D/investigation/scripts/epi_empirical_analysis_fun.R")

# load data
load("/Users/jpowell/repo/eQTL-2D/investigation/data/investigation_data.RData")

# set the data dir
setwd("/Users/jpowell/repo/eQTL-2D/investigation/data/output/")

# Read in 30 sig list
sig30 <- read.csv("/Users/jpowell/repo/eQTL-2D/investigation/data/sig_list.csv", header=T)

# list files
lf <- list.files()

# Read in the permutation results
pemp <- read.table("/Users/jpowell/repo/eQTL-2D/investigation/data/pemp_out.txt", header=T)


##################################################################
##################################################################
##################################################################
# summary of the clean and total snps

gs <- genome_sum.fun(lf)


##################################################################
##################################################################
##################################################################
# add the filter info and y/n in 30 sig

gs <- filter_add.fun(sig, gs)


##################################################################
##################################################################
##################################################################
# Determine the summary of the output

lambda <- summarize.fun(lf)
gs <- cbind(gs, lambda[,5:11])


##################################################################
##################################################################
##################################################################
# Genomic Control check
GC_out <- gc_check.fun(lf)

png(filename="~/repo/eQTL-2D/investigation/docs/figures/lambda.png")
plot(as.numeric(as.matrix(GC_out$lambdaC)), as.numeric(as.matrix(GC_out$lambdaF)), pch=16,
	xlab="Lambda GC - chisq", ylab="Lambda F")
dev.off()

# Merge GC_out and pemp data
l1 <- paste(pemp$probename, pemp$snp1, pemp$snp2, sep="_")
l2 <- paste(test$probename, test$snp1, test$snp2, sep="_")
index <- which(l2 %in% l1)
GC_out500 <- GC_out[index,]

# Lambda Chi-sq
png(filename="~/repo/eQTL-2D/investigation/docs/figures/perm_vs_lambdaC.png")
plot(as.numeric(as.matrix(pemp$neglog10pemp)), -log10(as.numeric(as.matrix(GC_out500$PlamC))),
	xlab="permutation", ylab="GWAS - lambda Chisq", main="", pch=16)
dev.off()


# Lambda F4
png(filename="~/repo/eQTL-2D/investigation/docs/figures/perm_vs_lambdaF.png")
plot(as.numeric(as.matrix(pemp$neglog10pemp)), -log10(as.numeric(as.matrix(GC_out500$PlamF))),
	xlab="permutation", ylab="GWAS - lambda F (4df)", main="", pch=16)
dev.off()



##################################################################
##################################################################
##################################################################
# Calculate the additive eQTL effects for each pair

add_test <- add_cal.fun(sig, bsgs, geno, bim)


##################################################################
##################################################################
##################################################################
# Make figures

png(filename="~/repo/eQTL-2D/investigation/docs/figures/lambda_npairs.png")
plot(as.numeric(as.matrix(gs$lambda)), as.numeric(as.matrix(gs$nthreshold)), 
		col=as.numeric(as.matrix(gs$filter)),
		pch=16,
		xlab="lambda", ylab="n snp pairs > 4.48x10-6")
dev.off()



# Make figures for the different filters
f1 <- which(gs$filter==1)
f2 <- which(gs$filter==2)

# F1
png(filename="~/repo/eQTL-2D/investigation/docs/figures/lambdaF1.png")
hist(as.numeric(as.matrix(gs$lambda[f1])), breaks=25, 
	xlab="lambda", main="Filter 1", col="lightgrey")
dev.off()


# F2
png(filename="~/repo/eQTL-2D/investigation/docs/figures/lambdaF2.png")
hist(as.numeric(as.matrix(gs$lambda[f2])), breaks=25, 
	xlab="lambda", main="Filter 2", col="lightgrey")
dev.off()



load(lf[3])
index <- which(output$nclass==9 & output$minclass > 5 & output$LD < 0.1)
foo <- output[index,]
png(filename="~/repo/eQTL-2D/investigation/docs/figures/lambda_107.png")
hist(foo$P, main="ILMN_1651886 - lambda = 1.07", xlab="p-value", col="lightgrey")
dev.off()


load(lf[17])
index <- which(output$nclass==9 & output$minclass > 5 & output$LD < 0.1)
foo <- output[index,]
png(filename="~/repo/eQTL-2D/investigation/docs/figures/lambda_305.png")
hist(foo$P, main="ILMN_1660549 - lambda = 3.05", xlab="p-value", col="lightgrey")
dev.off()

png(filename="~/repo/eQTL-2D/investigation/docs/figures/empPval.png")
hist(as.numeric(as.matrix(gs$P_emp)), breaks=20,
	main="", xlab="-log10 p-values", col="lightgrey")
dev.off()

png(filename="~/repo/eQTL-2D/investigation/docs/figures/empF.png")
hist(as.numeric(as.matrix(gs$F_emp)), breaks=20,
	main="", xlab="Empirical F-value", col="lightgrey")
dev.off()


# 1. For each probe-SNP pair (where SNP = fixed), calculate an empirical 5% type-I error rate, as you suggested yourself this morning. This will be useful to assess replication results. So count the number of observed test statistics greater than the 95th percentile of an F-distribution with [4,N-5] degrees of freedom. I'm not sure what the df in the denominator are, something like N-4 or N-5. 
png(filename="~/repo/eQTL-2D/investigation/docs/figures/Npairs_above_empthres.png")
hist(as.numeric(as.matrix(gs$N_F_empNtests)), breaks=25, 
	xlab="N pairs above emp 0.05 threshold", col="lightgrey", main="")
dev.off()

png(filename="~/repo/eQTL-2D/investigation/docs/figures/type1.png")
hist(as.numeric(as.matrix(gs$Type1)), breaks=35,
	xlab="Empirical type 1 error rate", col="lightgrey", main="")
dev.off()

# 2. Given results from 1.
# a) display the empirical type-I error rate for the probe-SNP-SNP trios that are in the Nature article table (the top 30).
gs30 <- type1_30.fun(gs, sig30)
xtable(gs30[,c(1:4,10,13,19, 23)])


# 3. For each probe-SNP pair, list 
# a) the largest observed F-statistic from the empirical results
# b) list the corresponding p-value (= 1/n I think, where n = number of SNPs that passed the filters)
# c) and list the F-statistic(s) from the probe-SNP-SNP trios that were among the 501

##################################################################
##################################################################
##################################################################
# 4. Calculate a mean lambda per probe, for those probes that were in the 501 selected trios. I am guessing that the lambdas for probes MBLN1 and TMEM149 are large.

multi_lambda <- multi_lambda.fun(gs,5)
xtable(multi_lambda)

# 5. For those probe-SNP-SNP trios that pass the empirical threshold (100 out of 400?), what proportion are cis and trans? I'm guessing that most will be cis.

##################################################################
##################################################################
##################################################################
# make histogram of the empirical p-values from permuation
png(filename="~/repo/eQTL-2D/investigation/docs/figures/pemp_hist.png")
hist(-log10(pemp$pemp), breaks=25, xlab="empirical (10,000,000) permutation p-value NegLog scale", main="", col="lightgray")
dev.off()

##################################################################
##################################################################
##################################################################

# Merge gs and pemp data
l1 <- paste(pemp$probename, pemp$snp1, pemp$snp2, sep="_")
l2 <- paste(gs$probename, gs$snp1, gs$snp2, sep="_")
index <- which(l2 %in% l1)
gs500 <- gs[index,]
n2 <- paste(gs500$probename, gs500$snp1, gs500$snp2, sep="_")

# P1
png(filename="~/repo/eQTL-2D/investigation/docs/figures/perm_vs_gwas.png")
plot(as.numeric(as.matrix(pemp$neglog10pemp)), as.numeric(as.matrix(gs500$P_emp)),
	xlab="permutation", ylab="GWAS", main="")
dev.off()

# P2
length(which(as.numeric(as.matrix(pemp$neglog10pemp)) > 5.34))
index <- which(as.numeric(as.matrix(pemp$neglog10pemp)) > 5.34)
tmp <- gs500[which(n2 %in% l1[index]),]
length(which(tmp$chr1!=tmp$chr2))



# P3
pemp30 <- type1_30.fun(pemp, sig30)
n1 <- paste(pemp30$probegene, pemp30$snp1, pemp30$snp2, sep="_")
n2 <- paste(gs500$gene, gs500$snp1, gs500$snp2, sep="_")

index <- which(as.numeric(as.matrix(pemp30$neglog10pemp)) > 5.34)
n3 <- n1[index]
index <- which(n2 %in% n3)
tmp <- gs500[index,]
xtable(tmp[,c(1:3,12:14,17,22)])


length(which(as.numeric(as.matrix(pemp30$neglog10pemp)) > 5.34))

#p4
f2 <- which(pemp$filter==2)
length(which(as.numeric(as.matrix(pemp30$neglog10pemp[f2])) > 5.34))









