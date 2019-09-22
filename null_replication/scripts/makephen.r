library(simulateGP)
library(dplyr)
library(data.table)

args <- commandArgs(T)

rawfile <- args[1]
famfile <- args[2]
varexp <- as.numeric(args[3])
polyfile <- args[4]
polyvar <- as.numeric(args[5])

print(rawfile)
print(famfile)
print(varexp)

gen <- fread(rawfile) %>% as_tibble()
names(gen)[ncol(gen)] <- "geno"

score <- fread(polyfile) %>% as_tibble()

stopifnot(all(score$IID == gen$IID))

gen$PHENOTYPE <- make_phen(c(sqrt(varexp), sqrt(polyvar)), cbind(gen$geno, score$SCORE))
cor(gen$PHENOTYPE, gen$geno, use="pair")^2
cor(gen$PHENOTYPE, score$SCORE, use="pair")^2

gen <- subset(gen, select=-c(geno))
write.table(gen, file=famfile, row=F, col=F, qu=F, sep=" ")



