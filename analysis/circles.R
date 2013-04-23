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


plotCircos <- function(gr, links, dot, anno=FALSE)
{
	if(anno)
	{
		a <- ggplot() + 
			layout_circle(gr, geom = "ideo", radius = 6, trackWidth = 1) +
			layout_circle(gr, geom = "text", radius = 6, trackWidth = 1, colour = "white", size=3, aes(label = seqnames)) +
			layout_circle(links, geom = "link", linked.to = "links2", radius = 3.5, trackwidth = 1) +
			layout_circle(dot, geom = "point", radius = 6, trackwidth = 1, colour="red", size=3, aes(y = score)) +
			layout_circle(dot, geom = "text", radius = 7, trackWidth = 1, colour = "black", size=2.5, aes(label = gene))
	} else {
		a <- ggplot() + 
			layout_circle(gr, geom = "ideo", radius = 6, trackWidth = 1) +
			layout_circle(links, geom = "link", linked.to = "links2", radius = 3.5, trackwidth = 1) +
			layout_circle(dot, geom = "point", radius = 6, trackwidth = 1, colour="red", size=3, aes(y = score)) +
			layout_circle(dot, geom = "text", radius = 7, trackWidth = 1, colour = "black", size=2.5, aes(label = gene))
	}
	return(a)
}


loadData <- function()
{
	load("~/repo/eQTL-2D/filtering/filtered_by_chr/sig_all_probes.RData")
	bim <- read.table("~/repo/eQTL-2D/data/clean_geno_final.bim", colClasses=c("character", "character", "numeric", "numeric", "character", "character"))
	sig$position1 <- bim$V4[sig$pos1]
	sig$position2 <- bim$V4[sig$pos2]
	sig$snp1 <- as.character(sig$snp1)
	sig$snp2 <- as.character(sig$snp2)
	sig$probename <- as.character(sig$probename)
	sig$probegene <- as.character(sig$probegene)
	return(sig)
}

sig <- loadData()
gr <- makeGr()
index <- table(sig$probename)
sig_mult <- subset(sig, probename %in% names(index)[index > 2])

a <- dlply(sig_mult, .(probeid), .progress="text", function(x)
{
	x <- mutate(x)
	links <- makeLinks(x, gr)
	dot <- makeDot(x, gr)
	a <- plotCircos(gr, links, dot)
	return(a)
})


multiplot(plotlist=a, cols=6)



