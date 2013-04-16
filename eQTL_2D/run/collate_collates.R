
n <- 20
start <- 1
end <- 2000
a <- seq(start, end, by=n)
b <- seq(n, end, by=n)

alldat <- data.frame()
for(i in 46:length(a))
{
    cat(i, "\n")
    file <- paste("res", a[i], "-", b[i], ".RData", sep="")
    load(file)
    alldat <- rbind(alldat, dat)
}





