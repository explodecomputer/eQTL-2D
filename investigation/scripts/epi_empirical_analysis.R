# epi_empirical_analysis.R

# analysis and interpretation of the output from epi_empirical.R
# joseph.powell@uq.edu.au

# read in the functions
source("/Users/jpowell/repo/eQTL-2D/investigation/scripts/epi_empirical_analysis_fun.R")

# load data
load("/Users/jpowell/repo/eQTL-2D/investigation/data/investigation_data.RData")

# set the data dir
setwd("/Users/jpowell/repo/eQTL-2D/investigation/data/output/")

# list files
lf <- list.files()


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
gs <- cbind(gs, lambda[,5:9])


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


png(filename="~/repo/eQTL-2D/investigation/docs/figures/lambda.png")
hist(as.numeric(as.matrix(gs$lambda)), breaks=25, 
	col=as.numeric(as.matrix(gs$filter)),
	xlab="lambda", main="")
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


hist(-log10(as.numeric(as.matrix(gs$P_emp))), breaks=20,
	main="", xlab="-log10 p-values", col="lightgrey")






