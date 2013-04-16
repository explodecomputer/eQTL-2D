# creating static data objects Ainv, fam, phedat
# remove IDs 8635101 and 8635102 (not present in genotype data)


setwd("wrayvisscher/eQTL_2D/data/")

# Create Ainv

library(pedigree)
ped <- read.table("ILMN_1720083.ped", header=T)
dim(ped)
remove <- c(8635101, 8635102)
remove <- which(ped[,2] %in% remove)
ped <- ped[-remove, ]
dim(ped)

fam <- data.frame(id=as.character(ped[,2]), dam=as.character(ped[,4]), sire=as.character(ped[,3]))
ord <- orderPed(fam)

#fam1 <- fam[order(ord), ]
#fam2 <- fam1[ord, ]
#head(fam)
#head(fam1)
#head(fam2)

fam <- fam[order(ord), ]
makeAinv(fam) # writes to Ainv.txt

Ai <- read.table("Ainv.txt")
Ainv <- matrix(0,nrow = nrow(ped),ncol = nrow(ped))
Ainv[as.matrix(Ai[,1:2])] <- Ai[,3]
dd <- diag(Ainv)
Ainv <- Ainv + t(Ainv)
diag(Ainv) <- dd
Ainv <- Ainv[ord,ord]

#y <- rnorm(nrow(ped))
#dim(y%*%Ainv2)


# Create phedat - matrix of phenotypes n x p

phedat <- read.csv("probe_signal_N10.csv", header=T) # remove rows 829 860
phedat <- phedat[, -c(remove, 829, 860)]
dim(phedat)

missing_count <- apply(phedat, 1, function(x) sum(is.na(x)))
#hist(missing_count)
#table(missing_count)

phedat <- phedat[missing_count == 0, ]
dim(phedat)
phedat <- t(matrix(unlist(phedat), nrow(phedat), ncol(phedat)))


# Create fam - first 6 columns of plink .ped file

fam <- ped[, c(1, 2, 3, 4, 5, 8)]


# save file

save(fam, Ainv, phedat, file="eqtl2d_objects.RData")

