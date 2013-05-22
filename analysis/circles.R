library(ggbio)
library(grid)
library(gridExtra)
library(plyr)
library(GenomicRanges)


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL)
{
	require(grid)

	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)

	numPlots = length(plots)

	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) 
	{
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
		ncol = cols, nrow = ceiling(numPlots/cols))
	}

	if (numPlots==1) 
	{
		print(plots[[1]])
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

		# Make each plot, in the correct location
		for (i in 1:numPlots)
		{
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
			layout.pos.col = matchidx$col))
		}
	}
}

ReadOrig <- function(filename, objname, setname)
{
	load(filename)
	sig <- get(objname)

	sig <- subset(sig, select=c(chr1, chr2, snp1, snp2, pos1, pos2, probename, probegene, probechr, pfull, pnest))
	sig$set <- setname
	sig$code <- with(sig, paste(probename, snp1, snp2))
	return(sig)
}

ReadRep <- function(filename, objname, setname)
{
	load(filename)
	sig <- get(objname)

	a <- sum(sig$replication_r^2 > 0.1)
	# cat(a, "pairs with high LD\n")

	b <- sum(sig$replication_nclass != 9)
	# cat(b, "pairs with < 9 classes\n")

	sig <- subset(sig, replication_r^2 < 0.1 & replication_nclass == 9, select=c(chr1, chr2, snp1, snp2, pos1, pos2, probename, probegene, probechr, replication_pfull, replication_pnest))

	names(sig) <- c("chr1", "chr2", "snp1", "snp2", "pos1", "pos2", "probename", "probegene", "probechr", "pfull", "pnest")
	sig$set <- setname
	sig$code <- with(sig, paste(probename, snp1, snp2))
	return(sig)
}

posData <- function(sig, bim)
{
	sig$position1 <- bim$V4[sig$pos1]
	sig$position2 <- bim$V4[sig$pos2]
	sig$snp1 <- as.character(sig$snp1)
	sig$snp2 <- as.character(sig$snp2)
	sig$probename <- as.character(sig$probename)
	sig$probegene <- as.character(sig$probegene)
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


makeGr <- function()
{
	data("hg19Ideogram", package = "biovizBase")
	chr.sub <- paste("chr", 1:22, sep = "")
	new.names <- as.character(1:22)
	names(new.names) <- paste("chr", new.names, sep = "")
	hg19Ideo <- hg19Ideogram
	hg19Ideo <- keepSeqlevels(hg19Ideogram, chr.sub)
	hg19Ideo <- renameSeqlevels(hg19Ideo, new.names)
	head(hg19Ideo)
	gr <- GRanges(
		seqnames = (1:22),
		ranges = hg19Ideo@ranges
	)
	seqlengths(gr) <- seqlengths(hg19Ideo)
	return(gr)
}


makeLinks <- function(x, gr)
{
	links1 <- GRanges(
		seqnames = x$chr1,
		IRanges(start = x$position1, width=1)
	)
	links2 <- GRanges(
		seqnames = x$chr2,
		IRanges(start = x$position2, width=1)
	)
	links1 <- c(gr, links1)
	links1 <- links1[-c(1:22)]
	links2 <- c(gr, links2)
	links2 <- links2[-c(1:22)]

	values(links1)$links2 <- links2


	return(links1)
}


makeDot <- function(x, gr)
{
	ran <- c(1, seqlengths(gr)[x$probechr[1]])
	names(ran) <- NULL
	dot <- GRanges(
		seqnames = x$probechr[1],
		IRanges(start = round(mean(ran)), width=1)
	)
	dot <- c(gr, dot)
	dot <- dot[-c(1:22)]
	values(dot)$score <- 0
	values(dot)$gene <- x$probegene[1]
	return(dot)
}


plotCircos <- function(gr, links, dot)
{
	a <- ggplot() + 
		layout_circle(gr, geom = "ideo", radius = 6, trackWidth = 1) +
		layout_circle(links, geom = "link", linked.to = "links2", radius = 3.5, trackwidth = 1) +
		layout_circle(dot, geom = "point", radius = 6, trackwidth = 1, colour="red", size=3, aes(y = score)) +
		layout_circle(dot, geom = "text", radius = 7, trackWidth = 1, colour = "black", size=2.5, aes(label = gene))

	return(a)
}


#=================================================================================================#
#=================================================================================================#


# Load data 

load("~/repo/eQTL-2D/analysis/replication_summary.RData")


# Circles 

gr <- makeGr()
index <- table(sig$probename)
sig_mult <- subset(sig, probename %in% names(index)[index > 12])


a <- dlply(sig_mult, .(probename), .progress="text", function(x)
{
	x <- mutate(x)
	links <- makeLinks(x, gr)
	dot <- makeDot(x, gr)
	a <- plotCircos(gr, links, dot)
	return(a)
})

plot(a[[1]])

multiplot(plotlist=a, cols=2)

