benjaminiHochberg <- function(pvalue, alpha)
{
	sorted.pvalue <- sort(pvalue[!is.na(pvalue)])
	n <- length(sorted.pvalue)
	j.alpha <- c(1:n) * (alpha / n)

	dif <- sorted.pvalue - j.alpha
	neg.dif <- dif[dif < 0]
	pos.dif <- neg.dif[length(neg.dif)]
	index <- dif == pos.dif
	p.cutoff <- sorted.pvalue[index]
	return(p.cutoff)
}


getType <- function(sig)
{
	type <- array("cis-cis", nrow(sig))
	type[sig$chr1 != sig$probechr & sig$chr2 == sig$probechr] <- "cis-trans"
	type[sig$chr1 == sig$probechr & sig$chr2 != sig$probechr] <- "cis-trans"
	type[sig$chr1 != sig$probechr & sig$chr2 != sig$probechr] <- "trans-trans"
	return(type)
}


#===========================================================================#
#===========================================================================#


load("~/eQTL-2D/analysis/interaction_list_meta_analysis.RData")
meta$type <- getType(meta)


sig <- subset(meta, filter != 3 & type != "cis-cis")
sig <- sig[order(sig$pnest_meta, decreasing=TRUE), ]
dim(sig)
names(sig)

p <- 10^-subset(sig, filter !=3)$pnest_meta
fdr <- benjaminiHochberg(p, 0.05)

a <- sig$pnest_meta > sig$upper

sig2 <- subset(sig, pnest_meta > -log10(fdr))
dim(sig2)

head(sig2)
subset(sig2)

sig2$type <- getType(sig2)
table(sig2$type, sig2$probegene)

