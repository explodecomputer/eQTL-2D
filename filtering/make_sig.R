# make sig

arguments <- commandArgs(T)
aggfile <- arguments[1]
T <- as.numeric(arguments[2])
alpha <- as.numeric(arguments[3])
outfile <- arguments[4]

load(aggfile)

sig <- subset(ld, pfull > T)
sig <- subset(sig, pnest > -log10(alpha/nrow(sig)))

SentinalSnp <- function(dat)
{
    dat$index <- 1:nrow(dat)
    dat$code <- with(dat, paste(chr1, chr2, probename))
    dat <- dat[order(dat$pnest), ]
    dat <- subset(dat, !duplicated(code))
    dat <- dat[order(dat$index), ]
    dat <- subset(dat, select=-c(index, code))
    return(dat)
}

sig <- SentinalSnp(sig)
dim(sig)

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

sig <- removeFactors(sig)

save(sig, file=outfile)

