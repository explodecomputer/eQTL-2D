
# For each of the significant interactions see if the two positions are within 1kb of any of the chromosome interactions

ciOverlap <- function(ci, sig, win)
{
	ci$index <- 1:nrow(ci)
	sig$int1 <- NA
	sig$int2 <- NA
	for(i in 1:nrow(sig))
	{
		cat(i, "of", nrow(sig), "\n")
		chr1 <- sig$chr1[i]
		chr2 <- sig$chr2[i]
		pos1 <- sig$pos1[i]
		pos2 <- sig$pos2[i]

		sa <- subset(ci, 
			loci1_chromosome == chr1 & 
			loci2_chromosome == chr2 &
			abs(loci1_position - pos1) <= win &
			abs(loci2_position - pos2) <= win
		)

		sb <- subset(ci, 
			loci1_chromosome == chr2 & 
			loci2_chromosome == chr1 &
			abs(loci1_position - pos2) <= win &
			abs(loci2_position - pos1) <= win
		)

		if(nrow(sa) > 0)
		{
			sig$int1[i] <- list(sa$index)
			print("Found!")
		}
		if(nrow(sb) > 0)
		{
			sig$int2[i] <- list(sb$index)
			print("Found!")
		}
	}
	return(sig)
}

createFake <- function(bim, n)
{
	a <- sample(1:nrow(bim), n, replace=FALSE)
	b <- sample(1:nrow(bim), n, replace=FALSE)

	dat <- data.frame(chr1=bim$V1[a], pos1=bim$V4[a], chr2=bim$V1[b], pos2=bim$V4[b])
	return(dat)
}


#=============================================================#
#=============================================================#


arguments <- commandArgs(T)

jid <- as.numeric(arguments[1])
n <- as.numeric(arguments[2])
win <- as.numeric(arguments[3])
cifile <- arguments[4]
bimfile <- arguments[5]
outroot <- arguments[6]

output <- paste(outroot, n, win, jid, sep="_")

#=============================================================#
#=============================================================#


ci <- read.csv(cifile, header=T)
bim <- read.table(bimfile, colClasses=c("numeric", "character", "numeric", "numeric", "character", "character"))

head(ci)
head(bim)

fake <- createFake(bim, n)
head(fake)

b <- ciOverlap(ci, fake, win)

write.table(sum(!is.na(b$int1) | !is.na(b$int2)), file=output, row=F, col=F, qu=F)
