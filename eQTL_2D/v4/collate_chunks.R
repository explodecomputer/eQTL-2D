rootname <- commandArgs(T)[1]
start <- as.numeric(commandArgs(T)[2])
end <- as.numeric(commandArgs(T)[3])
output <- commandArgs(T)[4]

alldat <- data.frame()
subdat <- data.frame()

for(i in start:end) {

    fn <- paste(rootname, i, ".RData", sep="")
    load(fn)
    alldat <- rbind(alldat, subset(dat, !is.na(chr1)))
    subdat <- rbind(subdat, subset(dat, propA <= 0.6))

}

save(alldat, file=paste(output, start, ".RData", sep=""))
save(subdat, file=paste(output, start, "_sub.RData", sep=""))


