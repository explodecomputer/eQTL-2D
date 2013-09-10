# zcat arichu_info.gz | awk '{if($4 > 0.05 && $4 < 0.95){print $0}}' | gzip > arichu_common_info.gz

info <- read.table("arichu_common_info.gz", header=F, colClasses="character")
info$V3 <- as.numeric(info$V3)
info$V4 <- as.numeric(info$V4)
info$V5 <- as.numeric(info$V5)
info$V6 <- as.numeric(info$V6)
info$V7 <- as.numeric(info$V7)
info$V8 <- as.numeric(info$V8)
info$V9 <- as.numeric(info$V9)
info$V10 <- as.numeric(info$V10)

a <- grep(":D", info$V2, value=TRUE)
b <- grep(":I", info$V2, value=TRUE)
length(a)
length(b)

info <- subset(info, ! V2 %in% c(a, b))

save(info, file="arichu_common_info.RData")
allsnps <- info$V2
save(allsnps, file="arichu_common_snps.RData")

gz1 <- gzfile("arichu_common_snps.txt.gz", "w")
write.table(allsnps, gz1, row=F, col=F, qu=F)
close(gz1)

