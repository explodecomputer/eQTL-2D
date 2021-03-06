rtdr <- commandArgs(T)[1]
output <- commandArgs(T)[2]
nom <- paste(rtdr, "filtered", 1:7400, ".RData", sep="")

l <- list()
j <- 1
for(i in 1:length(nom))
{
	if(file.exists(nom[i])){
		cat(i, "\n")
		load(nom[i])
		l[[j]] <- res
		j <- j+1
	}
}

library(plyr)
ld <- rbind.fill(l)

save(ld, file=output)
q()


head(ld)
hist(ld$pnest)hist(ld$propA)
hist(ld$pfull)
hist(ld$pint)
table(ld$pnest > -log10(0.05/nrow(ld)))
table(ld$pnest > -log10(0.05/4168), ld$pfull > 16.5)
table(table(ld$probename))

ldo <- ld[with(ld, order(probeid, snp1, pnest)), ]
ldo <- subset(ldo, !duplicated(paste(snp1, probename)))

ldo <- ld[with(ldo, order(probeid, snp2, pnest)), ]
ldo <- subset(ldo, !duplicated(paste(snp2, probename)))

table(table(ldo$probename))

hist(ldo$propA)

plot(propG ~ probehsq, data=ldo)
R --no-save --args 1 ../v4/results/result ../v4/scratch/resphen ../data/residuals.RData ../data/geno.RData 15.38 0.1 5 out < filter_raw_output.R
