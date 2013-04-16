
allsnps <- read.table("GWAS.SNPlist", header=F)
mapfile <- read.table("../../../info/GWAS.map", header=F)
mapfile$V4 <- 0

dim(mapfile)
snps <- as.character(allsnps[,1])
length(snps)

mapfile <- mapfile[mapfile[,2] %in% snps, ]
dim(mapfile)


for(i in 1:22)
{
    print(i)
    write.table(subset(mapfile, V1 == i), file=paste("GWAS_all.map", i, sep=""), row.names=F, col.names=F, quote=F)
}




