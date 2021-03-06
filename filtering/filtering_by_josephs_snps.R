load("filtered_by_chr/aggregated.RData")
ls()
dim(ld)
head(ld)
hist(ld$pint)
table(ld$pint > -log10(0.05 / nrow(ld))
)
a <- table(ld$pint > -log10(0.05 / nrow(ld)))
a <- subset(ld, pint > -log10(0.05 / nrow(ld)))
dim(a)
head(a)
table(a$pfull > 16.2)
table(table(a$probeid))
hist(a$pfull)
a <- subset(ld, pnest > -log10(0.05 / nrow(ld)))
dim(a)
head(a)
table(table(a$probeid))
table(a$probeid)
subset(a, probeid==4743)
hist(a$pnest)
with(a, table(duplicated(paste(chr1, chr2, probeid))))
table(a$pfull > 16.2)
b <- a[order(a$pnest), ]
head(b)
b <- a[order(a$pnest, decre=T), ]
b <- a[order(a$pnest, decreasing=T), ]
head(b)
hist(b$propA)
ls()
head(a)
head(b)
plot(pnest ~ propA, data=b)
ls()
head(a)
dim(b)
table(b$pfull > 16.2)
ls()
a <- read.table("Full_Merlin_out_NLP_8.txt", hea=T)
dim(a)
head(a)
a <- read.table("Full_Merlin_out_NLP_8.txt", hea=F)
dim(a)
head(a)
head(b)
ad <- paste(a$V2, a$V6)
dim(ad(
dim(ad)
length(ad)
head(ad)
ep1 <- paste(b$snp1, b$probename)
ep2 <- paste(b$snp2, b$probename)
ep <- unique(c(ep1, ep2))
length(ep)
length(ad)
table(ep %in% ad)
table(ad %in% ep)
epbig <- subset(ep, ad %in% ep)
length(epbig)
head(epbig)
table(is.na(ep))
table(is.na(epbig))
epbig <- subset(ep, ep %in% ad)
table(is.na(epbig))
epbig
ls()
head(b)
b$ep1 <- paste(b$snp1, b$probename)
b$ep2 <- paste(b$snp2, b$probename)
head(b)
head(subset(b, b$ep1 %in% epbig | b$ep2 %in% epbig))
a <- (subset(b, b$ep1 %in% epbig | b$ep2 %in% epbig))
head(a)
ls()
head(a)
head(ad)
head(b)
head(ad)
joseph <- read.table("Full_Merlin_out_NLP_8.txt")
head(joseph)
joseph$code <- paste(joseph$V2, joseph$V6)
j2 <- subset(joseph, code %in% epbig)
dim(j2)
head(j2)
hist(-log10(j2$V11))
hist(-log10(joseph$V11))
dev.new()
hist(-log10(j2$V11))
table(-log10(j2$V11) > 10)
table(-log10(joseph$V11) > 10)
epnew <- ep[! ep %in% ad]
length(epnew)
j3 <- subset(joseph, code %in% epnew)
dim(j3)
head(epnew)
head(epnew)
head(b)
head(epnew)
subset(b, ep1 %in% epnew & ep2 %in% epnew)
newepi <- subset(b, ep1 %in% epnew & ep2 %in% epnew)
dim(newepi)
newepi
dir()
getwd()
savehistory("filtering_by_josephs_snps.R")
