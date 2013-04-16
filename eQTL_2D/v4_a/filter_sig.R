args <- commandArgs(T)
aggfile <- args[1]
outfile <- args[2]

load(aggfile)

sig <- subset(aggfile, chr1 <= 22 & chr2 <= 22)
sig <- subset(aggfile, pnest > -log10(0.05/nrow(sig)))

dim(sig)
length(unique(c(sig$snp1, sig$snp2)))
length(unique(sig$probename))

save(sig, file=outfile)
