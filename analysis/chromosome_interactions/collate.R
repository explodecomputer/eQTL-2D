a <- array(NA, 10000)
win <- "1e+06"
n <- 549
for(i in 1:10000)
{
	a[i] <- scan(paste("results/out", n, win, i, sep="_"), what=numeric())
}
save(a, file=paste("collate_", n, "_",  win, ".RData",  sep=""))

