

# phenotype data is already corrected for age and sex. needs to be corrected for polygenic and covariates

#setwd("~/data/qimr/expression/")

famdat <- read.table("clean_geno_final.fam", header=F)
phendat <- read.table("probe_pheno_plink_final.txt", header=T)
covdat <- read.table("covariates_plink_final.txt", header=F)

phendat <- phendat[ phendat[, 2] %in% famdat[, 2], -c(1:2)]
index <- apply(phendat, 2, function(x) { ! any(x == -9) })
phendat <- phendat[, index]
dim(phendat)

save(famdat, phendat, file="clean_data_objects.RData")


# calculate A matrix

library(pedigree)

fam <- data.frame(id=as.character(famdat[,2]), dam=as.character(famdat[,4]), sire=as.character(famdat[,3]))
ord <- orderPed(fam)
fam <- fam[order(ord), ]

makeA(fam, which=rep(T,nrow(fam))) # writes to A.txt
A <- read.table("A.txt", header=F)

nInd <- nrow(fam)
Amat <- matrix(0,nrow = nInd,ncol = nInd)
Amat[as.matrix(A[,1:2])] <- A[,3]
dd <- diag(Amat)
Amat <- Amat + t(Amat)
diag(Amat) <- dd

Amat <- Amat[ord, ord]
Amat[lower.tri(Amat)] <- NA
A <- melt(Amat)
A <- A[!is.na(A[,3]), ]


A <- data.frame(A[,2], A[,1], 1000, A[,3])
write.table(A, "clean.grm", row=F, col=F, qu=F)
unlink("A.txt")
system("gzip clean.grm")

ids <- famdat[, 1:2]
write.table(ids, "clean.grm.id", row=F, col=F, qu=F)



# for each analysis 
## use existing covariates file
## create new fam file (for epiGPU)
## create new pheno file (for GCTA)
## run gcta --reml





