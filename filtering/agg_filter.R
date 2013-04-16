
library(plyr)
prune <- function(res, thresh)
{
    res$code <- with(res, paste(chr1, chr2))
    a <- ddply(res, .(code), function(x)
    {
	x <- mutate(x)
	m <- which(x$pint==max(x$pint))[1]
	return(x[m, ])
    })
    return(a)
}

a <- list()
j <- 1
for(i in 1:5400)
{
    if(file.exists(paste("filtered", i, ".RData", sep="")))
    {
	print(paste("filtered", i, ".RData", sep=""))
	load(paste("filtered", i, ".RData", sep=""))
	a[[j]] <- prune(res)
	j <- j+1
    }
}

set1 <- rbind.fill(a)
set2 <- subset(set1, pint > -log10(0.05/nrow(set1)))

dim(set1)
dim(set2)

save(set1, set2, file="agg_filtered.RData")

