library(ggplot2)

a <- read.table("../../investigation/data/total_analysis_data.txt", he=T)

pdf("figures/lambda.pdf")
hist(a$lambdaC, xlab="Genomic inflation", breaks=50, main=NULL, col="grey")
dev.off()
