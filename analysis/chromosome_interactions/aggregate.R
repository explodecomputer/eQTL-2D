rootname <- commandArgs(T)[1]
outname <- commandArgs(T)[2]

(cmd <- paste("cat ", rootname, "* > ", outname, sep=""))
system(cmd)

a <- scan(outname, what=numeric())
range(a)

save(a, file=paste(outname, ".RData", sep=""))

