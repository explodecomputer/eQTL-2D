library(ggbio)
library(grid)
library(gridExtra)
library(plyr)
library(GenomicRanges)
library(qgraph)

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
		IRanges(start = x$position1, width=1),
		seqinfo = seqinfo(gr)
	)
	links2 <- GRanges(
		seqnames = x$chr2,
		IRanges(start = x$position2, width=1),
		seqinfo = seqinfo(gr)
	)
	values(links1)$links2 <- links2
	values(links1)$col <- x$rep
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
		layout_circle(links, geom = "link", linked.to = "links2", radius = 3.5, trackwidth = 1, aes(colour=col)) +
		layout_circle(dot, geom = "point", radius = 6, trackwidth = 1, colour="red", size=3, aes(y = score)) +
		layout_circle(dot, geom = "text", radius = 7, trackWidth = 1, colour = "black", size=2.5, aes(label = gene)) +
		scale_colour_manual(values=c("light grey", "blue", "red"), drop = FALSE) +
		theme(legend.position="false")

	return(a)
}


linkColour <- function(sig)
{
	sig$rep <- 0
	sig$rep[sig$pnest_fehr > sig$upper_fehr | sig$pnest_egcut > sig$upper_egcut] <- 1
	sig$rep[sig$pnest_fehr > sig$upper_fehr & sig$pnest_egcut > sig$upper_egcut] <- 2
	sig$rep <- as.factor(sig$rep)
	return(sig)
}


multipleSnps <- function(sig)
{
	a <- subset(sig, select=c(snp1, snp2))
	b <- data.frame(t(apply(a, 1, sort)))
	names(b) <- c("snpa", "snpb")
	sig <- cbind(sig, b)
	sig$code2 <- with(sig, paste(snpa, snpb))
	tab <- table(sig$code2)
	nom <- names(tab)[tab>1]
	s <- subset(sig, code2 %in% nom)
	s <- s[order(s$code2), ]
}


#=================================================================================================#
#=================================================================================================#


# Load data 

load("~/repo/eQTL-2D/analysis/interaction_list_replication_summary.RData")
sig <- linkColour(sig)
sig <- subset(sig, filter != 3)

# Circles 

gr <- makeGr()
index <- table(sig$probename)
sig_mult <- subset(sig, probename %in% names(index)[index > 3])

thresh <- -log10(0.05 / 549)

sig_mult2 <- subset(sig, pnest_egcut > thresh | pnest_fehr > thresh)
sig_mult2$rep <- 1
sig_mult2$rep[sig_mult2$pnest_fehr > thresh & sig_mult2$pnest_egcut > thresh] <- 2
sig_mult2 <- subset(sig_mult2, !is.na(vc_fehr))
sig_mult2$rep <- as.factor(sig_mult2$rep)

# sig_mult <- subset(sig, probegene == "MBNL1")
# links <- makeLinks(sig_mult[1:13,], gr)

a <- dlply(sig_mult, .(probename), .progress="text", function(x)
{
	x <- mutate(x)
	links <- makeLinks(x, gr)
	dot <- makeDot(x, gr)
	a <- plotCircos(gr, links, dot)
	return(a)
})

b <- dlply(sig_mult2, .(probename), .progress="text", function(x)
{
	x <- mutate(x)
	links <- makeLinks(x, gr)
	dot <- makeDot(x, gr)
	a <- plotCircos(gr, links, dot)
	return(a)
})


pdf(file="~/repo/eQTL-2D/analysis/images/circles_replication2.pdf", width=25, height=20)
multiplot(plotlist=a, cols=6)
dev.off()


pdf(file="~/repo/eQTL-2D/analysis/images/circles_replication_bonf.pdf", width=25, height=20)
multiplot(plotlist=b, cols=6)
dev.off()



key <- subset(sig, !duplicated(probename), select=c(probename, probegene))
nom <- subset(key, duplicated(probegene))$probename

sig2 <- subset(sig, ! probename %in% nom)

# More circles
a <- subset(sig2, select=c(snp1, snp2))
names(a) <- c("from", "to")
a$thickness <- 1
a$col <- as.character(sig2$rep)
a$col[a$col == 0] <- "white"
a$col[a$col == 1] <- "black"
a$col[a$col == 2] <- "black"
a$col2 <- a$col
a$col2[sig2$rep == 1] <- "white"
a$col3 <- as.character(sig2$rep)
a$col3[a$col3 == 0] <- ""
a$col3[a$col3 == 1] <- "black"
a$col3[a$col3 == 2] <- "red"
a$thickness2 <- a$thickness
a$thickness2[a$col3 == "red"] <- 1.5
a$cistrans <- "blue"
a$cistrans[sig2$c]


pdf(file="~/repo/eQTL-2D/analysis/images/hairballs_all.pdf")
qgraph(a, curve=0.1, borders=FALSE, esize=1, color="grey", edge.color=a$thickness, labels=FALSE, vsize=0.2, arrows=FALSE)
dev.off()
pdf(file="~/repo/eQTL-2D/analysis/images/hairballs_rep.pdf")
qgraph(a, curve=0.1, borders=FALSE, esize=1, color="grey", edge.color=a$col, labels=FALSE, vsize=0.2, arrows=FALSE)
dev.off()
pdf(file="~/repo/eQTL-2D/analysis/images/hairballs_rep2.pdf")
qgraph(a, curve=0.1, borders=FALSE, esize=1, color="grey", edge.color=a$col2, labels=FALSE, vsize=0.2, arrows=FALSE)
dev.off()
pdf(file="~/repo/eQTL-2D/analysis/images/hairballs_allrep.pdf")
qgraph(a, curve=0.1, borders=FALSE, esize=a$thickness2, color="grey", edge.color=a$col3, labels=FALSE, vsize=0.2, arrows=FALSE)
dev.off()

pdf(file="~/repo/eQTL-2D/analysis/images/hairballs_all_reps.pdf", width=18, height=6)
par(mfrow=c(1,3))
qgraph(a, curve=0.1, borders=FALSE, esize=1, color="red", edge.color=a$thickness, labels=FALSE, vsize=0.3, arrows=FALSE)
qgraph(a, curve=0.1, borders=FALSE, esize=1, color="red", edge.color=a$col, labels=FALSE, vsize=0.2, arrows=FALSE)
qgraph(a, curve=0.1, borders=FALSE, esize=1, color="red", edge.color=a$col2, labels=FALSE, vsize=0.2, arrows=FALSE)
dev.off()


# Are there SNP pairs that affect more than one probe?
m <- multipleSnps(sig)

# in all cases this is because the different probes tag the same gene
