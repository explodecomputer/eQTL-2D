load("clean_geno_final.RData")
load("../analysis/interaction_list_meta_analysis.RData")

fam$index <- 1:nrow(fam)
fam_ind <- subset(fam, !duplicated(V1))
dim(fam)
dim(fam_ind)


snps <- with(subset(meta, filter!=3 & !duplicated(probegene)), unique(c(snp1)))

i <- fam_ind$index
j <- match(snps, bim$V2)

xmat[1:10,1:10]
xmat_ind <- xmat[i, j]
dim(xmat_ind)


r <- cor(xmat_ind, use="pair")^2
r[1:10,1:10]

l <- r[lower.tri(r, diag=F)]
length(l)
mean(l)
range(l)

1/length(i)

range(l)
hist(l)


