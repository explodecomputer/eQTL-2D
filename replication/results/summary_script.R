# Summary of the replication results
# Results sent by Harm-Jan Westra	

library(reshape2)
library(ggplot2)
library(VennDiagram)
library(xtable)


#'<brief desc>
#'
#'<full description>
#' @param filename <what param does>
#' @param  objname <what param does>
#' @param  setname <what param does>
#' @export
#' @keywords
#' @seealso
#' @return
#' @alias
#' @examples \dontrun{
#'
#'}
ReadOrig <- function(filename, objname, setname)
{
	load(filename)
	sig <- get(objname)

	sig <- subset(sig, select=c(chr1, chr2, snp1, snp2, probename, probegene, probechr, pfull, pnest))
	sig$set <- setname
	sig$code <- with(sig, paste(probename, snp1, snp2))
	return(sig)
}


#'<brief desc>
#'
#'<full description>
#' @param filename <what param does>
#' @param  objname <what param does>
#' @param  setname <what param does>
#' @export
#' @keywords
#' @seealso
#' @return
#' @alias
#' @examples \dontrun{
#'
#'}
ReadRep <- function(filename, objname, setname)
{
	load(filename)
	sig <- get(objname)

	a <- sum(sig$replication_r^2 > 0.1)
	cat(a, "pairs with high LD\n")

	b <- sum(sig$replication_nclass != 9)
	cat(b, "pairs with < 9 classes\n")

	sig <- subset(sig, replication_r^2 < 0.1 & replication_nclass == 9, select=c(chr1, chr2, snp1, snp2, probename, probegene, probechr, replication_pfull, replication_pnest))

	names(sig) <- c("chr1", "chr2", "snp1", "snp2", "probename", "probegene", "probechr", "pfull", "pnest")
	sig$set <- setname
	sig$code <- with(sig, paste(probename, snp1, snp2))
	return(sig)
}

confInt <- function(n, alpha)
{
	k <- c(1:n)
	upper <- -log10(qbeta(alpha/2,k,n+1-k))
	lower <- -log10(qbeta((1-alpha/2),k,n+1-k))
	expect <- -log10((k-0.5)/n)
	return(data.frame(expect, lower, upper))
}


qqDat <- function(sig, alpha)
{
	sig <- sig[order(sig$pnest, decreasing=T), ]
	ci <- confInt(nrow(sig), alpha)
	sig <- data.frame(sig, ci)
	return(sig)
}


panelCor <- function(x, y, digits=2, prefix="", cex.cor)
{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	r <- cor(x, y, use="pair")
	txt <- format(c(r, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor))
	{
		cex <- 0.8/strwidth(txt)
	}
	test <- cor.test(x,y)
	text(0.5, 0.5, txt, cex = cex)
}

uniteData <- function(..., pval)
{
	l <- list(...)

	# Get all pairs in common
	plist <- lapply(l, function(x) { return(x$code)})
	plist <- Reduce(intersect, plist)

	l <- lapply(l, function(x)
	{
		x <- subset(x, code %in% plist)
		x <- x[order(x$code), ]
	})

	a <- data.frame(code=l[[1]]$code)
	for(i in 1:length(l))
	{
		a <- data.frame(a, l[[i]][[pval]])
		names(a)[i+1] <- l[[i]]$set[1]
	}
	return(a)
}


plotCor <- function(ud, m=NULL)
{
	ud <- subset(ud, select=-c(code))
	ud <- as.matrix(ud)
	ud[is.infinite(ud)] <- NA
	a <- pairs(ud, lower.panel=panel.smooth, upper.panel=panelCor, main=m)
	return(a)
}


makeBonfTable <- function(fehr, egcut, marginal_list, top)
{
	# Choose the top 20 SNPs based on pnest from each set
	f <- fehr[order(fehr$pnest, decreasing=TRUE), ][1:top, ]
	e <- egcut[order(egcut$pnest, decreasing=TRUE), ][1:top, ]

	# Show the Bonferroni significant replications
	index <- f$pnest > -log10(0.05/nrow(fehr))
	f$set2 <- f$set
	f$set2[index] <- paste(f$set[index], "*", sep="")
	index <- e$pnest > -log10(0.05/nrow(fehr))
	e$set2 <- e$set
	e$set2[index] <- paste(e$set[index], "*", sep="")

	# Find the interactions common to both SNP sets
	index <- f$code %in% e$code
	nom <- f$code[index]
	f1 <- subset(f, code %in% nom)
	e1 <- subset(e, code %in% nom)
	f1 <- f1[order(f1$code), ]
	e1 <- e1[order(e1$code), ]
	f1$set2 <- paste(f1$set2, "+", e1$set2)
	fe <- rbind(f1, subset(e, ! code %in% nom), subset(f, ! code %in% nom))

	# Characterise the cis/trans types
	fe$type <- "cis-cis"
	fe$type[fe$chr1 == fe$probechr & fe$chr2 != fe$probechr] <- "cis-trans"
	fe$type[fe$chr2 == fe$probechr & fe$chr1 != fe$probechr] <- "trans-cis"
	fe$type[fe$chr2 != fe$probechr & fe$chr1 != fe$probechr] <- "trans-trans"


	# Characterise new SNPs vs known SNPs
	fe$code1 <- paste(fe$snp1, fe$probename)
	fe$code2 <- paste(fe$snp2, fe$probename)
	marginal_list$code <- paste(marginal_list$snp, marginal_list$probename)
	fe$margins <- "known-known"
	fe$margins[fe$code1 %in% marginal_list$code & ! fe$code2 %in% marginal_list$code] <- "known-new"
	fe$margins[! fe$code1 %in% marginal_list$code & fe$code2 %in% marginal_list$code] <- "new-known"
	fe$margins[! fe$code1 %in% marginal_list$code & ! fe$code2 %in% marginal_list$code] <- "new-new"


	# Tidy up the table
	fe$pg <- with(fe, paste(probegene, " (chr", probechr, ")", sep=""))
	fe$s1 <- with(fe, paste(snp1, " (chr", chr1, ")", sep=""))
	fe$s2 <- with(fe, paste(snp2, " (chr", chr2, ")", sep=""))

	fe <- subset(fe, select=c(pg, s1, s2, set2, type, margins))
	fe <- fe[order(fe$pg, fe$set2), ]
	names(fe) <- c("Probe gene", "SNP1", "SNP2", "Replication", "Type", "Previous associations")

	return(fe)
}


