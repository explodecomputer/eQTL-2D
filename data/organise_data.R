# Gib Hemani
# Read in pedigree data and create A matrix
# Read in phenotype data and create data object
# Test with ILMN_1653794 - should be ~ 0.5

library(GenABEL)
library(kinship)

fam <- read.table("clean_geno_final.fam", header=F)
dim(fam)
head(fam)

phendat <- read.table("probe_pheno_plink_final.txt", header=T)
phendat[1:10, 1:10]
dim(phendat)

covdat <- read.table("covariates_plink_final.txt", header=F)
head(covdat)
dim(covdat)

# Create A matrix

Amat <- kinship(fam[,2], fam[,3], fam[,4])
table(Amat)
Amat[1:10, 1:10]
Amat[upper.tri(Amat)] <- 1000000


# phendat - remove NA values, reorder

phendat <- phendat[phendat$IID %in% fam$V2, ]
dim(phendat)
all(phendat$IID == fam$V2)

missing <- apply(phendat, 2, function(x) any(x == -9))

phendat <- phendat[, ! missing]
dim(phendat)

id <- data.frame(IID=fam$V2)
phendat <- merge(id, phendat, by="IID")

# covariates

covdat <- covdat[ covdat$V2 %in% fam$V2, ]
head(covdat)
dim(covdat)

covdat <- merge(id, covdat, by.x="IID", by.y="V2")
head(covdat)
dim(covdat)

# testing

all(phendat$IID == fam$IID)
all(covdat$IID == fam$IID)


probe="ILMN_1653794"
ex <- data.frame(phen = phendat[, which(names(phendat) == probe)], fac1 = covdat$V3, fac2 = covdat$V4)
rownames(ex) <- id$IID

h <- polygenic(phen ~ fac1 + fac2, kin=Amat, ex)
h$esth2
# should equal 0.541


save(covdat, phendat, Amat, file="eqtl2d_objects.RData")


