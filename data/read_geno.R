# extract to 012 format using gcta then gunzip
# then tail -n +2 clean_geno_final.xmat | cut -d " " -f 3- > temp; mv temp clean_geno_final.xmat

xmatfile <- commandArgs(T)[1]
bimfile <- commandArgs(T)[2]
famfile <- commandArgs(T)[3]
outfile <- commandArgs(T)[4]

bim <- read.table(bimfile, header=F)
fam <- read.table(famfile, header=F)
nsnp <- nrow(bim)

xmat <- scan(xmatfile, what=integer())
nid <- length(xmat) / nsnp
xmat <- t(matrix(xmat, nsnp, nid))

save(xmat, bim, fam, file=outfile)




