

# phenotype data is already corrected for age and sex. needs to be corrected for polygenic and covariates

#setwd("~/data/qimr/expression/")

famdat <- read.table("clean_geno_final.fam", header=F)
phendat <- read.table("probe_pheno_plink_final.txt", header=T)
covdat <- read.table("covariates_plink_final.txt", header=F)

phendat <- phendat[ phendat[, 2] %in% famdat[, 2], -c(1:2)]


# calculate A matrix

library(pedigree)

fam <- data.frame(id=as.character(famdat[,2]), dam=as.character(famdat[,4]), sire=as.character(famdat[,3]))
ord <- orderPed(fam)
fam <- fam[order(ord), ]

makeA(fam, which=rep(T,nrow(fam))) # writes to A.txt
A <- read.table("A.txt", header=F)
A <- data.frame(A[,1], A[,2], 0, A[,3])
write.table(A, "clean.grm", row=F, col=F, qu=F)
unlink("A.txt")
system("gzip clean.grm")

ids <- famdat[ord, 1:2]
write.table(ids, "clean.grm.id", row=F, col=F, qu=F)


save(famdat, phendat, file="clean_data_objects.RData")

# for each analysis 
## use existing covariates file
## create new fam file (for epiGPU)
## create new pheno file (for GCTA)
## run gcta --reml





