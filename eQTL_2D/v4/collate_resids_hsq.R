# read in and store residuals and h2 estimates

rootname <- commandArgs(T)[1]
start <- as.numeric(commandArgs(T)[2])
end <- as.numeric(commandArgs(T)[3])
nid <- as.numeric(commandArgs(T)[4])
savefile <- commandArgs(T)[5]

n <- length(start:end)

resphen <- matrix(0, nid, n)
hsq <- array(0, n)

for(i in 1:n)
{
	cat(i,"\n")
	filename <- paste(rootname, i, sep="")
	hsqname <- paste(filename, ".h2", sep="")
	if(file.exists(filename))
	{
		resphen[, i] <- scan(filename)
	} else {
		cat("FILE MISSING: ", filename, "\n")
	}

	if(file.exists(hsqname))
	{
		hsq[i] <- scan(hsqname)
	} else {
		cat("FILE MISSING: ", hsqname, "\n")
	}
}

save(resphen, hsq, file=savefile)


