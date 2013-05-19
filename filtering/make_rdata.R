
read.egu <- function(rootname)
{
	if(!file.exists(paste(rootname, ".txt.gz", sep="")))
	{
		cat("Missing: ", rootname)
		return(data.frame(chr1=NA, chr2=NA, pos1=NA, pos2=NA,
			snp1=NA, snp2=NA, pfull=NA, pint=NA, df1=NA, df2=NA, complete=2))
	}

	complete <- system(paste("zgrep -q \"# 25 x 25 :\" ", rootname, ".txt.gz", sep=""))

	try(dat <- read.table(
		paste(rootname, ".txt.gz", sep=""), 
		skip=7, 
		header=T, 
		colClasses=c("numeric", "character", "numeric", "numeric", "character", "numeric", "numeric", "numeric", "numeric", "numeric")
	))
	cat(nrow(dat), "lines read\n")
	dat$pfull <- with(dat, -log10(pf(Fval, df1, df2, lower.tail=F)))
	dat$pint <- with(dat, -log10(pf(Fint, 4, df2, lower.tail=F)))

	if(nrow(dat) == 0)
	{
		cat(0)
		return(data.frame(chr1=NA, chr2=NA, pos1=NA, pos2=NA,
			snp1=NA, snp2=NA, pfull=NA, pint=NA, df1=NA, df2=NA, complete=complete))
	}

	dat <- with(dat, data.frame(chr1=Chr1, chr2=Chr2, pos1=SNP1, pos2=SNP2, 
		snp1=SNP1name, snp2=SNP2name, pfull, pint, df1, df2, complete=complete))

	a <- with(dat, paste(pos1, pos2))
	a <- duplicated(a)
	dat <- dat[!a, ]
	return(dat)
}

read.hsq <- function(rootname)
{
	if(!file.exists(paste(rootname, ".h2", sep="")))
	{
		cat("Missing hsq: ", rootname)
		return (NA)
	}
	a <- as.numeric(read.table(paste(rootname, ".h2", sep=""), header=F)[1,1])
	return(a)
}

read.phen <- function(rootname)
{
	if(!file.exists(rootname))
	{
		cat("Missing phen: ", rootname)
		return (NA)
	}
	a <- as.numeric(read.table(rootname, header=F)[,1])
	return(a)
}


# Get arguments
i <- as.numeric(commandArgs(T)[1])
rootres <- commandArgs(T)[2]
roothsq <- commandArgs(T)[3]
phenfile <- commandArgs(T)[4]
output <- commandArgs(T)[5]

# Load probe info
load(phenfile)
head(probeinfo)
probeinfo[i, ]

# Read in epiGPU file
res <- read.egu(paste(rootres, i, sep=""))
dim(res)
head(res)

# Add probe information
res$probeid <- i
res$probename <- probeinfo$PROBE_ID[i]
res$probechr <- probeinfo$CHROMOSOME_NEW[i]
res$probegene <- probeinfo$ILMN_GENE[i]

# Read hsq values
hsq <- read.hsq(paste(roothsq, i, sep=""))
res$probehsq <- hsq

head(res)

save(res, file=paste(output, i, ".RData", sep=""))
