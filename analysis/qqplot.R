library(reshape2)
library(ggplot2)
library(VennDiagram)
library(xtable)
library(grid)

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
	a <- subset(sig, filter==3)
	a <- a[order(a$pnest, decreasing=T), ]
	ci <- confInt(nrow(a), alpha)
	a <- data.frame(a, ci)

	b <- subset(sig, filter!=3)
	b <- b[order(b$pnest, decreasing=T), ]
	ci <- confInt(nrow(b), alpha)
	b <- data.frame(b, ci)

	return(rbind(b, a))
}


load("~/repo/eQTL-2D/analysis/interaction_list_replication_summary.RData")

fehr <- subset(sig_all, !is.na(pnest_fehr), select=c(filter, pnest_fehr))
egcut <- subset(sig_all, !is.na(pnest_egcut), select=c(filter, pnest_egcut))
names(fehr) <- names(egcut) <- c("filter", "pnest")
fehr <- qqDat(data.frame(fehr, set="Fehrmann"), 0.05)
egcut <- qqDat(data.frame(egcut, set="EGCUT"), 0.05)


dat1 <- rbind(fehr, egcut)
dat1$ex <- FALSE
dat1$ex[dat1$pnest > -log10(0.05/(nrow(dat1)/2))] <- TRUE
dat1$fake <- "Null"
dat1$fake[dat1$filter != 3] <- "Observed"

qqplot1 <- ggplot(subset(dat1, filter != 3)) +
	geom_ribbon(aes(x=expect, ymin=lower, ymax=upper), colour="white", alpha=0.5) +
	geom_abline(intercept=0, slope=1) +
	geom_point(aes(x=expect, y=pnest, colour=ex)) +
	geom_hline(yintercept=-log10(0.05/(nrow(dat1)/2)), linetype="dotted") +
	facet_grid(. ~ set) +
	labs(colour="Above 5% Bonferroni?", y="Observed", x="Expected") +
	theme(legend.position="none")
ggsave(qqplot1, file="~/repo/eQTL-2D/analysis/images/qqbonf.pdf", height=3.5, width=7)

# Show FDR correction
dat2 <- rbind(fehr[-c(1:20), ], egcut[-c(1:20), ])
dat2$ex <- FALSE
dat2$ex[dat2$pnest > dat2$upper] <- TRUE
dat2$fake <- "Null"
dat2$fake[dat2$filter != 3] <- "Observed"

qqplot2 <- ggplot(dat2) +
	geom_ribbon(aes(x=expect, ymin=lower, ymax=upper), colour="white", alpha=0.5) +
	geom_abline(intercept=0, slope=1) +
	geom_point(aes(x=expect, y=pnest, colour=ex), size=1.5) +
	facet_grid(set ~ fake) +
	labs(colour="Above FDR 5% CI?", y="Observed", x="Expected") +
	scale_colour_brewer(type="qual", palette=3) +
	theme(legend.position="none")
ggsave(qqplot2, file="~/repo/eQTL-2D/analysis/images/qqfdr.pdf")

