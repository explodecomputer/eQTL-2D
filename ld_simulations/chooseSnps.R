# R --no-save --args /clusterdata/uqgheman/ibimp/arichu/data/target/chr\*/ARICHU\*

ar <- commandArgs(T)

chipsnpsfile <- ar[1]
idlistfile <- ar[2]
plinkrt <- ar[3]
nqtl <- as.numeric(ar[4])

keepids <- scan(idlistfile, what="character")
chipsnps <- scan(chipsnpsfile, what="character")

getAllSnps <- function(plinkrt)
{
	require(plyr)
	bim <- list()
	for(i in 1:22)
	{
		cat(i, "\n")
		filename <- gsub("\\*", i, plinkrt)
		bim[[i]] <- data.frame(chr = i, rsid = read.table(filename, colClasses="character")$V2)
	}
	bim <- rbind.fill(bim)
	return(bim)
}




