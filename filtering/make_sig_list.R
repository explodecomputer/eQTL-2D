
removeFactors <- function(x)
{
    nom <- names(x)
    l <- nom[lapply(x, class)=="factor"]
    for(i in 1:length(l))
    {
	x[[l[i]]] <- as.character(x[[l[i]]])
    }
    return(x)
}


load("nonmarginal_interaction_list.RData")
sig1 <- removeFactors(sig)
load("old/sig.RData")
sig <- removeFactors(sig)


sig$code <- with(sig, paste(probename, snp1, snp2))
sig1$code <- with(sig1, paste(probename, snp1, snp2))

table(sig1$code %in% sig$code)
table(sig$code %in% sig1$code)

subset(sig, ! code %in% sig1$code)

sig <- rbind(sig, sig1)
sig <- subset(sig, !duplicated(code))
dim(sig)

table(table(sig$probename))

signm <- sig
signm$filter <- 2

load("../../replication/run/interactions_list.RData")
sig <- removeFactors(sig)
sig$filter <- 1

sig <- rbind(sig, subset(signm, select=-code))
sig$code <- with(sig, paste(probename, snp1, snp2))
table(duplicated(sig$code), sig$filter)

sig <- subset(sig, !duplicated(code))


# Get random list

load("~/repo/eQTL-2D/data/residuals_all.RData")
bim <- read.table("~/repo/eQTL-2D/data/clean_geno_final.bim", colClasses="character")

sig2 <- sig
sig2$filter <- 3
bim <- subset(bim, V1 %in% 1:22)
index1 <- sample(1:nrow(bim), nrow(sig2))
index2 <- sample(1:nrow(bim), nrow(sig2))
table(index1 %in% index2)

sig2$chr1 <- as.numeric(bim$V1[index1])
sig2$snp1 <- bim$V2[index1]
sig2$pos1 <- index1
sig2$chr2 <- as.numeric(bim$V1[index2])
sig2$snp2 <- bim$V2[index2]
sig2$pos2 <- index2
sig2$pfull <- sig2$pint <- sig2$df1 <- sig2$df2 <- sig2$minclasssize <- sig2$propG <- sig2$propA <- NA

sig2$code <- with(sig2, paste(probename, snp1, snp2))
head(sig2)
head(sig)

sig <- rbind(sig, sig2)
sig <- subset(sig, select=-code)

save(sig, file="~/repo/eQTL-2D/replication/run/interactions_list2.RData")



