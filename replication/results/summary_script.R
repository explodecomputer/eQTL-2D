# Summary of the replication results
# Results sent by Harm-Jan Westra	

library(reshape2)
library(ggplot2)
library(VennDiagram)



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



# Load Data

bsgs <- ReadOrig("~/repo/eQTL-2D/replication/run/interactions_list.RData", "sig", "BSGS")
egcut <- ReadRep("~/repo/eQTL-2D/replication/results/EGCUT_replication.RData", "newsig", "EGCUT")
fehr <- ReadRep("~/repo/eQTL-2D/replication/results/FehrmannHT12v3_replication.RData", "newsig", "Ferhmann")


# Q-Q confidence intervals
egcut <- qqDat(egcut, 0.05)
fehr <- qqDat(fehr, 0.05)


# Bonferroni corrected significant replication
egcut_b <- subset(egcut, pnest > -log10(0.05 / nrow(egcut)))
fehr_b <- subset(fehr, pnest > -log10(0.05 / nrow(fehr)))


# FDR corrected significant replication
egcut_f <- subset(egcut, pnest > upper)
fehr_f <- subset(fehr, pnest > upper)


# Q-Q plot

# Show bonferroni correction
dat1 <- rbind(fehr, egcut)
dat1$ex <- FALSE
dat1$ex[dat1$pnest > -log10(0.05/(nrow(dat1)/2))] <- TRUE

qqplot1 <- ggplot(dat1) +
	geom_ribbon(aes(x=expect, ymin=lower, ymax=upper), colour="white", alpha=0.5) +
	geom_abline(intercept=0, slope=1) +
	geom_point(aes(x=expect, y=pnest, colour=ex)) +
	geom_hline(yintercept=-log10(0.05/(nrow(dat1)/2)), linetype="dotted") +
	facet_grid(. ~ set) +
	labs(colour="Above 5% Bonferroni?", y="Observed", x="Expected")


# Show FDR correction
dat2 <- rbind(fehr[-c(1:20), ], egcut[-c(1:20), ])
dat2$ex <- FALSE
dat2$ex[dat2$pnest > dat2$upper] <- TRUE

qqplot2 <- ggplot(dat2) +
	geom_ribbon(aes(x=expect, ymin=lower, ymax=upper), colour="white", alpha=0.5) +
	geom_abline(intercept=0, slope=1) +
	geom_point(aes(x=expect, y=pnest, colour=ex), size=1.5) +
	facet_grid(. ~ set) +
	labs(colour="Above FDR 5% CI?", y="Observed", x="Expected")



# Characterise overlap
# Significant in only Fehr, Significant in only EGCUT, Significant in both

# Make venn diagram


bonf.venn <- draw.triple.venn(
	area1 = nrow(bsgs),
	area2 = nrow(fehr_b),
	area3 = nrow(egcut_b),
	n12   = sum(fehr_b$code %in% bsgs$code),
	n23  = sum(egcut_b$code %in% fehr_b$code),
	n13  = sum(egcut_b$code %in% bsgs$code),
	n123 = sum(egcut_b$code %in% fehr_b$code),
	category=c("BSGS", "Fehrmann", "EGCUT"),
	fill = c("blue", "red", "green"),
	alpha=c(0.1, 0.1, 0.1),
	lty = "blank",
	cex = 2,
	cat.cex = 2,
	cat.col = c("blue", "red", "green"),
	ind=FALSE
)

fdr.venn <- draw.triple.venn(
	area1 = nrow(bsgs),
	area2 = nrow(fehr_f),
	area3 = nrow(egcut_f),
	n12   = sum(fehr_f$code %in% bsgs$code),
	n23  = sum(egcut_f$code %in% fehr_f$code),
	n13  = sum(egcut_f$code %in% bsgs$code),
	n123 = sum(egcut_f$code %in% fehr_f$code),
	category=c("BSGS", "Fehrmann", "EGCUT"),
	fill = c("blue", "red", "green"),
	alpha=c(0.1, 0.1, 0.1),
	lty = "blank",
	cex = 2,
	cat.cex = 2,
	cat.col = c("blue", "red", "green"),
	ind=FALSE
)


# Plot correlations of pfull and pnest

pfull_all <- uniteData(bsgs, egcut, fehr, pval="pfull")
pfull_f <- uniteData(bsgs, egcut_f, fehr_f, pval="pfull")
pnest_all <- uniteData(bsgs, egcut, fehr, pval="pnest")
pnest_f <- uniteData(bsgs, egcut_f, fehr_f, pval="pnest")

plotCor(pfull_all)
dev.new()
plotCor(pfull_f)

plotCor(pnest_all)
dev.new()
plotCor(pnest_f)



# List of top interactions

makeBonfTable <- function(fehr, egcut, marginal_list)
{
	f <- fehr[order(fehr$pnest, decreasing=TRUE), ][1:20, ]
	e <- egcut[order(egcut$pnest, decreasing=TRUE), ][1:20, ]

	fe <- rbind(f,e)
	fe$set[duplicated(fe$code)] <- "a"
	fe <- fe[order(fe$set), ]
	fe <- subset(fe, !duplicated(code))
	fe$set[fe$set=="a"] <- "Fehrmann + EGCUT"
	fe <- fe[order(fe$probegene, fe$set), ]
	fe$pg <- with(fe, paste(probegene, " (chr", probechr, ")", sep=""))
	fe$s1 <- with(fe, paste(snp1, " (chr", chr1, ")", sep=""))
	fe$s2 <- with(fe, paste(snp2, " (chr", chr2, ")", sep=""))
	fe$type <- "cis-cis"
	fe$type[fe$chr1 == fe$probechr & fe$chr2 != fe$probechr] <- "cis-trans"
	fe$type[fe$chr2 == fe$probechr & fe$chr1 != fe$probechr] <- "trans-cis"
	fe$type[fe$chr2 != fe$probechr & fe$chr1 != fe$probechr] <- "trans-trans"



	fe <- subset(fe, select=c(pg, s1, s2, set, type))

}



# Which SNPs are novel / have previously known marginal effects


# Which circles replicate












