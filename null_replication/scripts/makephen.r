library(simulateGP)
library(dplyr)

args <- commandArgs(T)

rawfile <- "../data/discsen.raw"
famfile <- "../data/disc.fam2"
rawfile <- args[1]
famfile <- args[2]
varexp <- as.numeric(args[3])

gen <- read_delim(rawfile, " ", col_names=TRUE)
names(gen)[ncol(gen)] <- "geno"

gen$PHENOTYPE <- make_phen(sqrt(varexp), gen$geno)
cor(gen$PHENOTYPE, gen$geno, use="pair")^2
gen <- subset(gen, select=-c(geno))
write.table(gen, file=famfile, row=F, col=F, qu=F, sep=" ")

