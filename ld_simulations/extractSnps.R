# Usage:

# R --no-save --args snplist.txt /ibscratch/wrayvisscher/imputation/arichu/data/imputed/chr\*/arichu_1kg_p1v3_\* outputfile < extractSnps.R

# - The first argument is the SNP list (just a list of rs SNP ids and nothing else)
# - The second argument specifies the root name for the genotype data (in binary plink format). Replace the chromosome number with "\*"
# - The third argument is the output file root name that you want to save everything to.


library(snpStats)
library(plyr)

findChromosomes <- function(snp, plinkrt)
{
	cmd <- paste("grep -P -m 1 '", snp, "\t' ", plinkrt, ".bim | head -n 1 | cut -d \":\" -f 2 | cut -f 1", sep="")
	cat(snp)
	chr <- system(cmd, intern=TRUE)
	cat(" ", chr, "\n")
	return(chr)
}


findChromosomesAll <- function(snplist, plinkrt)
{
	chr <- snplist
	for(i in 1:length(snplist))
	{
		a <- findChromosomes(snplist[i], plinkrt)
		chr[i] <- ifelse(is.null(a), NA, a)
	}
	dat <- data.frame(snp = snplist, chr = as.numeric(chr))
	dat <- subset(dat, !is.na(chr))
	dat$snp <- as.character(dat$snp)
	dat <- dat[order(dat$chr), ]
	return(dat)
}


findSnps <- function(snplist, bimfile)
{
	bim <- scan(bimfile, what="character")
	snps <- bim[seq(2, length(bim), 6)]
	keepsnps <- snps[snps %in% snplist]
	return(keepsnps)
}

findSnpsAll <- function(snplist, plinkrt)
{
	require(plyr)
	l <- list()
	for(i in 1:22)
	{
		cat(i, ": ")
		bimfile <- paste(gsub("\\*", i, plinkrt), ".bim", sep="")
		l[[i]] <- data.frame(snp = findSnps(snplist, bimfile), chr = i)
		cat(nrow(l[[i]]), "\n")
	}
	l <- rbind.fill(l)
	return(l)
}

extractSnps <- function(snpnames, plinkrt)
{
	require(snpStats)
	rawdata <- read.plink(bed=plinkrt, select.snps=snpnames)
	return(rawdata)
}

extractSnpsAll <- function(snpdat, plinkrt)
{
	chr <- unique(snpdat$chr)
	l <- length(chr)

	cat(chr[1], ":", snpdat$snp[snpdat$chr == chr[1]], "\n")
	dat <- extractSnps(snpdat$snp[snpdat$chr == chr[1]], gsub("\\*", chr[1], plinkrt))
	if(l == 1) return(dat)
	for(i in 2:l)
	{
		cat(chr[i], ":", snpdat$snp[snpdat$chr == chr[i]], "\n")
		dat1 <- extractSnps(snpdat$snp[snpdat$chr == chr[i]], gsub("\\*", chr[i], plinkrt))
		dat$map <- rbind(dat$map, dat1$map)
		dat$genotypes <- cbind(dat$genotypes, dat1$genotypes)
	}
	return(dat)
}

extractInfo <- function(snp, plinkrt)
{
	cmd <- paste("zgrep -P '", snp, " ' ", plinkrt, " | head -n 1 | cut -d \":\" -f 2", sep="")
	info <- system(cmd, intern=TRUE)
	print(info)
	return(info)
}

extractInfoAll <- function(snpdat, plinkrt)
{
	l <- nrow(snpdat)
	cmd <- paste("zcat ", gsub("\\*", "1", plinkrt), "_info.txt.gz | head -n 1", sep="")
	print(cmd)
	dat <- system(cmd, intern=TRUE)
	for(i in 1:l)
	{
		nom <- paste(gsub("\\*", snpdat$chr[i], plinkrt), "_info.txt.gz", sep="")
		snp <- snpdat$snp[i]
		cat(i, ":", snp, "\n")
		dat <- c(dat, extractInfo(snp, nom))
	}
	return(dat)	
}


ar <- commandArgs(T)
snplistfile <- ar[1]
plinkrt <- ar[2]
output <- ar[3]


snplist <- scan(snplistfile, what="character")
# snpdat <- findChromosomesAll(snplist, plinkrt)
# info <- extractInfoAll(snpdat, plinkrt)
snpdat <- findSnpsAll(snplist, plinkrt)
save(snpdat, file="temp.RData")
dat <- extractSnpsAll(snpdat, plinkrt)


# write.table(info, file=paste(output, "_info.txt", sep=""), row=F, col=F, qu=F)
write.plink(file.base = output, 
	snps = dat$genotypes, 
	pedigree = dat$fam$pedigree,
	id = dat$fam$member,
	father = dat$fam$father,
	mother = dat$fam$mother,
	sex = dat$fam$sex,
	phenotype = dat$fam$affected,
	chromosome = dat$map$chromosome, 
	genetic.distance = dat$map$cM, 
	position = dat$map$position, 
	allele.1 = dat$map$allele.1, 
	allele.2 = dat$map$allele.2
)

