# Gib Hemani
# Replicate epistatic signals

library(GenABEL)

args       <- commandArgs(T)
plink      <- args[1]
plinkfile  <- args[2]
signalfile <- args[3]
ourfile    <- args[4]
# covfile    <- args[5]
outfile    <- "replication.RData"





#################
# Organise data #
#################



# Check for files

stopifnot(file.exists(plink))
stopifnot(file.exists(plinkfile))
stopifnot(file.exists(signalfile))
# stopifnot(file.exists(covfile))
stopifnot(file.exists(ourfile))


# Read in data
signal <- read.table(signalfile, header=T)
rownames(signal) <- signal$IID
signal <- subset(signal, select=-c(IID, FID))

# covs   <- read.table(covfile, header=T)
load(ourfile)


# List of SNPs to extract from plink
snplist <- with(sig, unique(c(as.character(snp1), as.character(snp2))
write.table(snplist, file="snplist.txt", row=F, col=F, qu=F)


# Extract to plink format
cmd <- paste(plink, " --bfile ", plinkfile, " --extract snplist.txt --recode --out ", plinkfile, sep="")
stat <- system(cmd, intern=FALSE)
stopifnot(stat==0)

# Convert plink to 0/1/2 format

plink_to_012 <- function(pedfile, mapfile)
{
	pformat <- read.table(pedfile, header=FALSE, colClasses="character")
	map <- read.table(mapfile, header=FALSE, colClasses=c("character", "character", "numeric", "numeric"))
	ids <- pformat[, 1:6]
	nid <- nrow(ids)
	pformat <- pformat[, -c(1:6)]
	index <- seq(1, ncol(pformat), 2)
	geno <- matrix(0, nid, length(index))
	for(i in 1:length(index))
	{
		snp <- index[, c(index[i], index[i]+1)]
		x <- array(NA, nid)
		snp[snp == "0"] <- NA
		alleles <- names(sort(table(snp)))
		i0 <- snp[,1] == alleles[1] & snp[,2] == alleles[1]
		i2 <- snp[,1] == alleles[2] & snp[,2] == alleles[2]
		i1 <- (snp[,1] == alleles[1] & snp[,2] == alleles[2]) | (snp[,1] == alleles[2] & snp[,2] == alleles[1])
		x[i0] <- 0
		x[i1] <- 1
		x[i2] <- 2
		geno[, i] <- x
	}
	colnames(geno) <- map$V2
	rownames(geno) <- ids$V2
	return(geno)
}

geno <- plink_to_012(paste(plinkfile, ".ped", sep=""), paste(plinkfile, ".map", sep=""))
ids <- rownames(geno)
snps <- colnames(geno)
probes <- colnames(signal)[-c(1,2)] # First two columns should be FID and IID

# Check IDs are all unique
stopifnot(!any(duplicated(ids)))

# Check IDs are present in cov and signal files
# common_ids <- ids[ids %in% covs$IID & ids %in% signal$IID]
common_ids <- ids[ids %in% rownames(signal)]
stopifnot(length(common_ids) > 500)

geno <- geno[rownames(geno) %in% common_ids, ]
signal <- signal[rownames(signal) %in% common_ids, ]
# covs <- covs[covs$IID %in% common_ids, ]

stopifnot(all(rownames(geno) %in% rownames(signal)))

# reorder rows to match each other
index <- match(rownames(signal), rownames(geno))
geno <- geno[index, ]

stopifnot(all(rownames(geno) == rownames(signal)))




##################################
# Ready to run statistical tests #
##################################


sig$replication_pfull <- NA
sig$replication_pint <- NA
sig$replication_pnest <- NA
sig$replication_p1 <- NA
sig$replication_p2 <- NA
sig$replication_r <- NA
sig$replication_nclass <- NA
sig$replication_nid <- NA


for(i in 1:nrow(sig))
{
	cat(i, "of", nrow(sig), "\n")
	snp1 <- geno[colnames(geno) == sig$snp1]
	snp2 <- geno[colnames(geno) == sig$snp2]
	probe <- signal[colnames(signal) == sig$probename]

	sig$replication_p1[i] <- mean(snp1) / 2
	sig$replication_p2[i] <- mean(snp2) / 2
	sig$replication_r[i] <- cor(snp1, snp2, use="pair")
	sig$replication_nclass[i] <- length(table(snp1 + 3*snp2))
	sig$replication_nid[i] <- table(!is.na(snp1) & !is.na(snp2))
	fullmod <- lm(probe ~ as.factor(snp1)*as.factor(snp2))
	intmod <- lm(probe ~ as.factor(snp1):as.factor(snp2))
	temp <- summary(fullmod)$fstatistic
	sig$replication_pfull[i] <- -log10(pf(temp[1], temp[2], temp[3], low=FALSE))

	temp <- summary(fullmod)$fstatistic
	sig$replication_pfull[i] <- -log10(pf(temp[1], temp[2], temp[3], low=FALSE))
}



















