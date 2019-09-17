library(simulateGP)
library(dplyr)
library(data.table)

args <- commandArgs(T)

rawfile <- args[1]
famfile <- args[2]
varexp <- as.numeric(args[3])

print(rawfile)
print(famfile)
print(varexp)

gen <- fread(rawfile) %>% as_tibble()
names(gen)[ncol(gen)] <- "geno"

gen$PHENOTYPE <- make_phen(sqrt(varexp), gen$geno)
cor(gen$PHENOTYPE, gen$geno, use="pair")^2
gen <- subset(gen, select=-c(geno))
write.table(gen, file=famfile, row=F, col=F, qu=F, sep=" ")

