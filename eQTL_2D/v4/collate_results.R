
# results files are massive. this script will trim and collate them into RData files

datstat <- function(info, xmat, resphen)
{
    cat(" ... ",nrow(info), "rows to be analysed\n")
    a <- array(0, nrow(info))
    g <- array(0, nrow(info))
    for(i in 1:nrow(info))
    {
	if(i %% (nrow(info) / 100) == 0) cat(i/(nrow(info)/100)," ")
	l <- info[i,]
	if(is.na(l$chr1)) {
	    g[i] <- NA
	    a[i] <- NA
	    next
	}
	mod <- anova(lm(resphen[,l$probeid] ~ xmat[,l$pos1] + xmat[,l$pos2] + as.factor(xmat[,l$pos1]) : as.factor(xmat[,l$pos2])))
	g[i] <- sum(mod$Sum[1:3]) / mod$Sum[4]
	a[i] <- sum(mod$Sum[1:2]) / sum(mod$Sum[1:3])

    }
    return(data.frame(g,a))
}

datstat2 <- function(info, xmat, resphen)
{
    cat(nrow(info), "rows to be analysed\n")
    a <- array(0, nrow(info))
    g <- array(0, nrow(info))
    for(i in 1:nrow(info))
    {
        if(i %% (nrow(info) / 100) == 0) cat(i/(nrow(info)/100)," ")
        l <- info[i,]
        if(is.na(l$chr1)) {
            g[i] <- NA
            a[i] <- NA
            next
        }

        mod <- anova(lm(resphen[,l$probeid] ~ as.factor(xmat[,l$pos1]) * as.factor(xmat[,l$pos2])))
	print(mod)

	print(tapply(resphen[,l$probeid], list(xmat[,l$pos1], xmat[,l$pos2]), mean))
	print(table(xmat[,l$pos1], xmat[,l$pos2]))

    }
}




read.egu <- function(rootname, threshold)
{
    if(!file.exists(paste(rootname, ".txt", sep="")))
    {
	if(file.exists(paste(rootname, ".txt.gz", sep="")))
	{
	    system(paste("gunzip", paste(rootname, ".txt.gz", sep="")))
	    read.egu(rootname, threshold, test)
	}
	cat("Missing: ", rootname)
	return(data.frame(chr1=NA, chr2=NA, pos1=NA, pos2=NA,
            snp1=NA, snp2=NA, pfull=NA, pint=NA, df1=NA, df2=NA, complete=2))
    }

    complete <- system(paste("grep -q \"# 25 x 25 :\" ", rootname, ".txt", sep=""))

    dat <- read.table(paste(rootname, ".txt", sep=""), skip=7, header=T)
    dat <- subset(dat, df1 > 5)
    dat$pfull <- with(dat, -log10(pf(Fval, df1, df2, lower.tail=F)))
    dat$pint <- with(dat, -log10(pf(Fint, 4, df2, lower.tail=F)))
    dat <- subset(dat, pfull >= threshold | pint >= threshold)
    dat <- subset(dat, (Chr1 != Chr2) | (abs(SNP1 - SNP2) > 5))

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
   # cat(nrow(dat))
    return(dat)
}

#a <- read.egu("results/result1", 12)

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

######################

rootres <- commandArgs(T)[1]
phenfile <- commandArgs(T)[2]
genofile <- commandArgs(T)[3]
roothsq <- commandArgs(T)[4]
first <- commandArgs(T)[5]
last <- commandArgs(T)[6]
threshold <- as.numeric(commandArgs(T)[7])
eqtlobj <- commandArgs(T)[8]
output <- commandArgs(T)[9]


# probe | probe number | hsq | complete | chr1 | chr2 | pos1 | pos2 | snp1 | snp2 | Pfull | Pint | df1 | df2 
# if probe has no values then single row with NA 

load(phenfile)
load(genofile)
load(eqtlobj)

dat <- data.frame()
nid <- nrow(phendat)

for(i in first:last)
{
    cat(i, "... ")
    res <- read.egu(paste(rootres, i, sep=""), threshold)
    hsq <- read.hsq(paste(roothsq, i, sep=""))
    
    probeid <- names(phendat)[i+2]
    res$probe <- probeid
    res$probeid <- i
    res$hsq <- hsq
    temp <- datstat(res, xmat, resphen)
    res$propG <- temp$g
    res$propA <- temp$a

    dat <- rbind(dat, res)
    cat("\n")
}

save(dat, file=output)


